/**
 * MoleculeEditor - Main molecule editor component
 * 
 * Phase 2: Centralized editor state and operations
 * 
 * This component manages:
 * - Molecule state (using Molecule class from Phase 1)
 * - Drawing operations (add atom, add bond, delete, move)
 * - Selection state
 * - Hover state
 * - Tool mode
 */

import React, { useState, useCallback, useRef, useEffect } from 'react'
import { Molecule, validateMolecule, canCreateBond, generate2DLayout } from '@/lib/molecule'
import type { EditorTool, Atom, Bond, ValidationResult } from '@/lib/molecule'
import { CanvasLayer } from './CanvasLayer'
import { PointerManager, KeyboardManager } from '@/lib/molecule/input'
import { HistoryManager } from '@/lib/molecule/history'
import {
  AddAtomCommand,
  RemoveAtomCommand,
  AddBondCommand,
  RemoveBondCommand,
  MoveAtomCommand,
  ClearMoleculeCommand,
} from '@/lib/molecule/history'
import { nanoid } from 'nanoid'

interface MoleculeEditorProps {
  width: number
  height: number
  initialMolecule?: Molecule
  tool?: EditorTool
  onMoleculeChange?: (molecule: Molecule) => void
  onAtomSelect?: (atomId: string | null) => void
  onBondSelect?: (bondId: string | null) => void
}

export function MoleculeEditor({
  width,
  height,
  initialMolecule,
  tool = 'select',
  onMoleculeChange,
  onAtomSelect,
  onBondSelect,
}: MoleculeEditorProps) {
  const [molecule, setMolecule] = useState<Molecule>(initialMolecule || new Molecule())
  const [selectedAtomId, setSelectedAtomId] = useState<string | null>(null)
  const [selectedBondId, setSelectedBondId] = useState<string | null>(null)
  const [hoveredAtomId, setHoveredAtomId] = useState<string | null>(null)
  const [hoveredBondId, setHoveredBondId] = useState<string | null>(null)
  const [elementToAdd, setElementToAdd] = useState<string>('C')
  const [bondOrder, setBondOrder] = useState<number>(1)
  const [scale, setScale] = useState(1)
  const [offsetX, setOffsetX] = useState(0)
  const [offsetY, setOffsetY] = useState(0)
  
  // Input managers
  const pointerManagerRef = useRef<PointerManager>(new PointerManager())
  const keyboardManagerRef = useRef<KeyboardManager>(new KeyboardManager())
  
  // History manager
  const historyManagerRef = useRef<HistoryManager>(new HistoryManager())
  
  // Bond creation state (drag from atom to atom)
  const [bondStartAtomId, setBondStartAtomId] = useState<string | null>(null)
  
  // History state (for UI feedback)
  const [canUndo, setCanUndo] = useState(false)
  const [canRedo, setCanRedo] = useState(false)
  
  // Validation state
  const [validationResult, setValidationResult] = useState<ValidationResult | null>(null)

  // Update molecule and notify parent
  const updateMolecule = useCallback((newMolecule: Molecule) => {
    setMolecule(newMolecule)
    
    // Validate molecule
    const validation = validateMolecule(newMolecule)
    setValidationResult(validation)
    
    onMoleculeChange?.(newMolecule)
  }, [onMoleculeChange])
  
  // Validate molecule on mount and when it changes
  useEffect(() => {
    if (molecule) {
      const validation = validateMolecule(molecule)
      setValidationResult(validation)
    }
  }, [molecule])

  // Add atom at position
  const addAtom = useCallback((x: number, y: number, element: string = elementToAdd) => {
    // Convert screen coordinates to world coordinates
    const worldX = (x - width / 2 - offsetX) / scale
    const worldY = (y - height / 2 - offsetY) / scale
    const worldZ = 0

    const atomId = nanoid()
    const command = new AddAtomCommand({ id: atomId, element, position: [worldX, worldY, worldZ] })
    
    const updated = historyManagerRef.current.execute(molecule, command)
    updateMolecule(updated)
    setCanUndo(historyManagerRef.current.canUndo())
    setCanRedo(historyManagerRef.current.canRedo())
    
    return atomId
  }, [molecule, width, height, scale, offsetX, offsetY, elementToAdd, updateMolecule])

  // Add bond between two atoms (with validation)
  const addBond = useCallback((atom1Id: string, atom2Id: string, order: number = bondOrder) => {
    if (atom1Id === atom2Id) return null

    // Validate bond creation
    const validation = canCreateBond(molecule, atom1Id, atom2Id, order)
    if (!validation.valid) {
      // Could emit error event here
      console.warn('Cannot create bond:', validation.error?.message)
      return null
    }

    const bondId = nanoid()
    const command = new AddBondCommand(bondId, atom1Id, atom2Id, order)
    
    const updated = historyManagerRef.current.execute(molecule, command)
    updateMolecule(updated)
    setCanUndo(historyManagerRef.current.canUndo())
    setCanRedo(historyManagerRef.current.canRedo())
    
    return bondId
  }, [molecule, bondOrder, updateMolecule])

  // Delete atom (and all its bonds)
  const deleteAtom = useCallback((atomId: string) => {
    const command = new RemoveAtomCommand(atomId)
    
    const updated = historyManagerRef.current.execute(molecule, command)
    updateMolecule(updated)
    setCanUndo(historyManagerRef.current.canUndo())
    setCanRedo(historyManagerRef.current.canRedo())
    
    if (selectedAtomId === atomId) {
      setSelectedAtomId(null)
      onAtomSelect?.(null)
    }
  }, [molecule, selectedAtomId, updateMolecule, onAtomSelect])

  // Delete bond
  const deleteBond = useCallback((bondId: string) => {
    const command = new RemoveBondCommand(bondId)
    
    const updated = historyManagerRef.current.execute(molecule, command)
    updateMolecule(updated)
    setCanUndo(historyManagerRef.current.canUndo())
    setCanRedo(historyManagerRef.current.canRedo())
    
    if (selectedBondId === bondId) {
      setSelectedBondId(null)
      onBondSelect?.(null)
    }
  }, [molecule, selectedBondId, updateMolecule, onBondSelect])

  // Move atom (with history tracking)
  const moveAtom = useCallback((atomId: string, x: number, y: number, isDrag: boolean = false) => {
    const worldX = (x - width / 2 - offsetX) / scale
    const worldY = (y - height / 2 - offsetY) / scale
    const atom = molecule.getAtom(atomId)
    if (!atom) return

    // For drag operations, we batch moves into a single command
    // For now, create a command for each move (can be optimized later)
    const command = new MoveAtomCommand(atomId, [worldX, worldY, atom.position[2]])
    
    historyManagerRef.current.execute(molecule, command)
    updateMolecule(molecule)
    setCanUndo(historyManagerRef.current.canUndo())
    setCanRedo(historyManagerRef.current.canRedo())
  }, [molecule, width, height, scale, offsetX, offsetY, updateMolecule])

  // Handle atom click
  const handleAtomClick = useCallback((atomId: string, event: MouseEvent) => {
    if (tool === 'select') {
      setSelectedAtomId(atomId)
      setSelectedBondId(null)
      onAtomSelect?.(atomId)
      onBondSelect?.(null)
    } else if (tool === 'bond') {
      if (bondStartAtomId && bondStartAtomId !== atomId) {
        // Complete bond creation
        addBond(bondStartAtomId, atomId)
        setBondStartAtomId(null)
      } else {
        // Start bond creation
        setBondStartAtomId(atomId)
        setSelectedAtomId(atomId)
        onAtomSelect?.(atomId)
      }
    } else if (tool === 'delete') {
      deleteAtom(atomId)
    } else if (tool === 'add-atom') {
      // Place atom at click position (handled by canvas click)
    }
  }, [tool, bondStartAtomId, addBond, deleteAtom, onAtomSelect, onBondSelect])

  // Handle bond click
  const handleBondClick = useCallback((bondId: string, event: MouseEvent) => {
    if (tool === 'select') {
      setSelectedBondId(bondId)
      setSelectedAtomId(null)
      onBondSelect?.(bondId)
      onAtomSelect?.(null)
    } else if (tool === 'delete') {
      deleteBond(bondId)
    }
  }, [tool, deleteBond, onBondSelect, onAtomSelect])

  // Handle canvas click (for adding atoms)
  const handleCanvasClick = useCallback((x: number, y: number) => {
    if (tool === 'add-atom') {
      addAtom(x, y)
    }
  }, [tool, addAtom])

  // Clear molecule
  const clear = useCallback(() => {
    const command = new ClearMoleculeCommand()
    
    historyManagerRef.current.execute(molecule, command)
    updateMolecule(molecule)
    setCanUndo(historyManagerRef.current.canUndo())
    setCanRedo(historyManagerRef.current.canRedo())
    
    setSelectedAtomId(null)
    setSelectedBondId(null)
    setBondStartAtomId(null)
    onAtomSelect?.(null)
    onBondSelect?.(null)
  }, [molecule, updateMolecule, onAtomSelect, onBondSelect])

  // Undo
  const undo = useCallback(() => {
    const result = historyManagerRef.current.undo(molecule)
    if (result) {
      updateMolecule(result)
      setCanUndo(historyManagerRef.current.canUndo())
      setCanRedo(historyManagerRef.current.canRedo())
    }
  }, [molecule, updateMolecule])

  // Redo
  const redo = useCallback(() => {
    const result = historyManagerRef.current.redo(molecule)
    if (result) {
      updateMolecule(result)
      setCanUndo(historyManagerRef.current.canUndo())
      setCanRedo(historyManagerRef.current.canRedo())
    }
  }, [molecule, updateMolecule])

  // Setup keyboard shortcuts
  useEffect(() => {
    const keyboard = keyboardManagerRef.current

    // Escape - cancel current operation
    const unsubCancel = keyboard.on((event) => {
      if (event.key === 'Escape') {
        setBondStartAtomId(null)
        setSelectedAtomId(null)
        setSelectedBondId(null)
        onAtomSelect?.(null)
        onBondSelect?.(null)
        pointerManagerRef.current.cancel()
      }
    })

    // Delete/Backspace - delete selected
    const unsubDelete = keyboard.on((event) => {
      if (event.key === 'Delete' || event.key === 'Backspace') {
        if (selectedAtomId) {
          deleteAtom(selectedAtomId)
        } else if (selectedBondId) {
          deleteBond(selectedBondId)
        }
      }
    })

    // Undo/Redo shortcuts
    keyboard.register({
      key: 'z',
      ctrl: true,
      meta: true, // Support Cmd on Mac
      handler: (event) => {
        if (event.shiftKey) {
          redo()
        } else {
          undo()
        }
      },
      description: 'Undo (Ctrl+Z) / Redo (Ctrl+Shift+Z)',
    })

    // Setup keyboard event listener
    const handleKeyDown = (e: KeyboardEvent) => {
      keyboard.processEvent(e as any)
    }

    window.addEventListener('keydown', handleKeyDown)
    return () => {
      window.removeEventListener('keydown', handleKeyDown)
      unsubCancel()
      unsubDelete()
    }
  }, [selectedAtomId, selectedBondId, deleteAtom, deleteBond, undo, redo, onAtomSelect, onBondSelect])

  return (
    <div className="molecule-editor" style={{ width, height, position: 'relative' }}>
      <CanvasLayer
        molecule={molecule}
        selectedAtomId={selectedAtomId}
        selectedBondId={selectedBondId}
        hoveredAtomId={hoveredAtomId}
        hoveredBondId={hoveredBondId}
        width={width}
        height={height}
        tool={tool}
        scale={scale}
        offsetX={offsetX}
        offsetY={offsetY}
        onAtomClick={handleAtomClick}
        onBondClick={handleBondClick}
        onAtomHover={setHoveredAtomId}
        onBondHover={setHoveredBondId}
        pointerManager={pointerManagerRef.current}
        bondStartAtomId={bondStartAtomId}
      />
    </div>
  )
}

// Export editor operations for external use
export type MoleculeEditorRef = {
  addAtom: (x: number, y: number, element?: string) => string
  addBond: (atom1Id: string, atom2Id: string, order?: number) => string | null
  deleteAtom: (atomId: string) => void
  deleteBond: (bondId: string) => void
  moveAtom: (atomId: string, x: number, y: number) => void
  clear: () => void
  undo: () => void
  redo: () => void
  canUndo: () => boolean
  canRedo: () => boolean
  getMolecule: () => Molecule
  autoLayout: () => Promise<void>
}

// Export undo/redo state for UI
export interface MoleculeEditorState {
  canUndo: boolean
  canRedo: boolean
  undoDescription: string | null
  redoDescription: string | null
}

