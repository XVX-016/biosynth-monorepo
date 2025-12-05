/**
 * EditorContext - Master context for molecule editor
 * 
 * Phase 10: Lab Page UI Rebuild
 * 
 * Unifies:
 * - molecule state
 * - selection state
 * - current tool
 * - predictions
 * - history
 */

import React, { createContext, useContext, useState, useCallback, useRef, useEffect, useMemo } from 'react'
import { Molecule, predictionService } from '@/lib/molecule'
import type { EditorTool, ValidationResult, PredictionResult } from '@/lib/molecule'
import { HistoryManager } from '@/lib/molecule/history'
import { validateMolecule } from '@/lib/molecule/validation/Validator'
import { setupAutosave, loadFromLocalStorage } from '@/lib/molecule/storage/autosave'
import type { AttentionData } from '@/lib/molecule/prediction'
import { mapAttentionToMolecule } from '@/lib/molecule/attention'
import type { AttentionMap } from '@/lib/molecule/attention'

interface EditorContextValue {
  // Molecule state
  molecule: Molecule
  setMolecule: (molecule: Molecule) => void
  
  // Selection
  selectedAtomId: string | null
  selectedBondId: string | null
  setSelectedAtomId: (id: string | null) => void
  setSelectedBondId: (id: string | null) => void
  
  // Tool
  tool: EditorTool
  setTool: (tool: EditorTool) => void
  
  // Bond order (for bond tool)
  bondOrder: number
  setBondOrder: (order: number) => void
  
  // Validation
  validationResult: ValidationResult | null
  
  // Predictions
  predictions: PredictionResult | null
  predictionsLoading: boolean
  predictionsError: string | null
  
  // History
  canUndo: boolean
  canRedo: boolean
  undo: () => void
  redo: () => void
  
  // Stable SMILES (for predictions)
  stableSmiles: string | null

  // Attention overlay
  attentionData: AttentionData | null
  attentionMap: AttentionMap | null
  setAttentionData: (data: AttentionData | null) => void
  attentionOverlayEnabled: boolean
  setAttentionOverlayEnabled: (enabled: boolean) => void

  // Model selection
  selectedModelId: string
  setSelectedModelId: (modelId: string) => void
}

const EditorContext = createContext<EditorContextValue | null>(null)

export function useEditorContext() {
  const context = useContext(EditorContext)
  if (!context) {
    throw new Error('useEditorContext must be used within EditorProvider')
  }
  return context
}

interface EditorProviderProps {
  initialMolecule?: Molecule
  children: React.ReactNode
}

export function EditorProvider({ initialMolecule, children }: EditorProviderProps) {
  // Try to load from LocalStorage if no initial molecule provided
  const savedMolecule = initialMolecule || loadFromLocalStorage() || new Molecule()
  const [molecule, setMoleculeState] = useState<Molecule>(savedMolecule)
  const [selectedAtomId, setSelectedAtomId] = useState<string | null>(null)
  const [selectedBondId, setSelectedBondId] = useState<string | null>(null)
  const [tool, setTool] = useState<EditorTool>('select')
  const [bondOrder, setBondOrder] = useState<number>(1)
  const [validationResult, setValidationResult] = useState<ValidationResult | null>(null)
  const [predictions, setPredictions] = useState<PredictionResult | null>(null)
  const [predictionsLoading, setPredictionsLoading] = useState(false)
  const [predictionsError, setPredictionsError] = useState<string | null>(null)
  const [stableSmiles, setStableSmiles] = useState<string | null>(null)
  const [attentionDataState, setAttentionDataState] = useState<AttentionData | null>(null)
  const [attentionOverlayEnabled, setAttentionOverlayEnabled] = useState<boolean>(true)
  const [selectedModelId, setSelectedModelId] = useState<string>('propnet')
  
  const historyManagerRef = useRef(new HistoryManager())
  const predictionsTimeoutRef = useRef<NodeJS.Timeout | null>(null)

  // Update molecule and validate
  const setMolecule = useCallback((newMolecule: Molecule) => {
    setMoleculeState(newMolecule)
    
    // Validate
    const validation = validateMolecule(newMolecule)
    setValidationResult(validation)
    
    // Update history state
    setCanUndo(historyManagerRef.current.canUndo())
    setCanRedo(historyManagerRef.current.canRedo())
  }, [])

  const [canUndo, setCanUndo] = useState(false)
  const [canRedo, setCanRedo] = useState(false)

  // Undo/Redo
  const undo = useCallback(() => {
    const result = historyManagerRef.current.undo(molecule)
    if (result) {
      setMolecule(result)
    }
  }, [molecule, setMolecule])

  const redo = useCallback(() => {
    const result = historyManagerRef.current.redo(molecule)
    if (result) {
      setMolecule(result)
    }
  }, [molecule, setMolecule])

  const attentionMap = useMemo(() => {
    if (!attentionDataState || molecule.isEmpty()) {
      return null
    }
    try {
      return mapAttentionToMolecule(molecule, attentionDataState)
    } catch (error) {
      console.error('Failed to map attention data:', error)
      return null
    }
  }, [attentionDataState, molecule])

  const setAttentionData = useCallback((data: AttentionData | null) => {
    setAttentionDataState(data)
  }, [])

  // Debounced prediction updates
  useEffect(() => {
    if (molecule.isEmpty() || !validationResult?.valid) {
      setPredictions(null)
      setStableSmiles(null)
       setAttentionData(null)
      return
    }

    // Clear previous timeout
    if (predictionsTimeoutRef.current) {
      clearTimeout(predictionsTimeoutRef.current)
    }

    // Debounce predictions (300ms)
    predictionsTimeoutRef.current = setTimeout(async () => {
      setPredictionsLoading(true)
      setPredictionsError(null)

      try {
        // Generate SMILES for stable prediction
        const { toSMILES } = await import('@/lib/molecule/export')
        const smiles = await toSMILES(molecule, true)
        setStableSmiles(smiles)

        // Get predictions
        const result = await predictionService.predictWithAttention(molecule, selectedModelId)
        setPredictions(result)
        setAttentionData(result.attention || null)
      } catch (error) {
        setPredictionsError(error instanceof Error ? error.message : 'Prediction failed')
        setPredictions(null)
        setAttentionData(null)
      } finally {
        setPredictionsLoading(false)
      }
    }, 300)

    return () => {
      if (predictionsTimeoutRef.current) {
        clearTimeout(predictionsTimeoutRef.current)
      }
    }
  }, [molecule, validationResult, selectedModelId, setAttentionData])

  // Validate on molecule change
  useEffect(() => {
    const validation = validateMolecule(molecule)
    setValidationResult(validation)
  }, [molecule])

  // Setup autosave
  useEffect(() => {
    if (!molecule.isEmpty()) {
      const cleanup = setupAutosave(molecule)
      return cleanup
    }
  }, [molecule])

  const value: EditorContextValue = {
    molecule,
    setMolecule,
    selectedAtomId,
    selectedBondId,
    setSelectedAtomId,
    setSelectedBondId,
    tool,
    setTool,
    bondOrder,
    setBondOrder,
    validationResult,
    predictions,
    predictionsLoading,
    predictionsError,
    canUndo,
    canRedo,
    undo,
    redo,
    stableSmiles,
    attentionData: attentionDataState,
    attentionMap,
    setAttentionData,
    attentionOverlayEnabled,
    setAttentionOverlayEnabled,
    selectedModelId,
    setSelectedModelId,
  }

  return <EditorContext.Provider value={value}>{children}</EditorContext.Provider>
}

