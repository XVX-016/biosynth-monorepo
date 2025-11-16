import { create } from 'zustand'
import { MoleculeGraph } from '@biosynth/engine'
import type { Element } from '@biosynth/engine'
import { predict, generate, predictFast } from '../lib/api'

interface BackendPredictions {
  stability?: number
  toxicity?: number
  solubility?: number
  bioavailability?: number
  novelty?: number
}

export type Tool = 'select' | 'add-atom' | 'bond' | 'delete'

interface MoleculeState {
  // State
  currentMolecule: MoleculeGraph | null
  autoBond: boolean
  selectedAtomId: string | null
  selectedBondId: string | null
  loadingState: 'idle' | 'loading' | 'success' | 'error'
  backendPredictions: BackendPredictions | null
  error: string | null
  tool: Tool
  atomToAdd: Element | null
  currentBondOrder: number

  // Actions
  setMolecule: (molecule: MoleculeGraph | null) => void
  selectAtom: (id: string | null) => void
  selectBond: (id: string | null) => void
  updatePosition: (id: string, position: [number, number, number]) => void
  fetchPredictions: () => Promise<void>
  fetchPredictionsFast: () => Promise<void>
  generateMolecule: (prompt: string) => Promise<void>
  reset: () => void
  setTool: (tool: Tool) => void
  resetTool: () => void
  setAtomToAdd: (element: Element | null) => void
  setBondOrder: (order: number) => void
  setAutoBond: (value: boolean) => void
}

const initialState = {
  currentMolecule: null,
  autoBond: true,
  selectedAtomId: null,
  selectedBondId: null,
  loadingState: 'idle' as const,
  backendPredictions: null,
  error: null,
  tool: 'select' as Tool,
  atomToAdd: null as Element | null,
  currentBondOrder: 1,
}

export const useMoleculeStore = create<MoleculeState>((set, get) => ({
  ...initialState,

  setMolecule: (molecule) => {
    set({ currentMolecule: molecule, error: null })
  },

  selectAtom: (id) => {
    set({ selectedAtomId: id })
  },

  selectBond: (id) => {
    set({ selectedBondId: id })
  },

  updatePosition: (id, position) => {
    const { currentMolecule } = get()
    if (!currentMolecule) return

    const atom = currentMolecule.atoms.get(id)
    if (!atom) return

    currentMolecule.atoms.set(id, { ...atom, position })
    set({ currentMolecule: currentMolecule.clone() })
  },

  fetchPredictions: async () => {
    const { currentMolecule } = get()
    if (!currentMolecule) return

    set({ loadingState: 'loading', error: null })

    try {
      // TODO: Use MoleculeSerializer.toSMILES() when implemented
      const smiles = 'C' // Placeholder
      const result = await predict(smiles)
      
      set({
        loadingState: 'success',
        backendPredictions: result.properties,
        error: null,
      })
    } catch (error) {
      set({
        loadingState: 'error',
        error: error instanceof Error ? error.message : 'Failed to fetch predictions',
        backendPredictions: null,
      })
    }
  },

  fetchPredictionsFast: async () => {
    const { currentMolecule } = get()
    if (!currentMolecule) return

    set({ loadingState: 'loading', error: null })

    try {
      // TODO: Use MoleculeSerializer.toSMILES() when implemented
      const smiles = 'C' // Placeholder
      const result = await predictFast(smiles)
      
      set({
        loadingState: 'success',
        backendPredictions: result.properties,
        error: null,
      })
    } catch (error) {
      set({
        loadingState: 'error',
        error: error instanceof Error ? error.message : 'Failed to fetch predictions',
        backendPredictions: null,
      })
    }
  },

  generateMolecule: async (prompt: string) => {
    set({ loadingState: 'loading', error: null })

    try {
      const result = await generate(prompt)
      // TODO: Use MoleculeSerializer.fromSMILES() when implemented
      // For now, create a placeholder molecule
      const molecule = new MoleculeGraph()
      molecule.addAtom({ element: 'C', position: [0, 0, 0] })
      
      set({
        currentMolecule: molecule,
        loadingState: 'success',
        error: null,
      })
    } catch (error) {
      set({
        loadingState: 'error',
        error: error instanceof Error ? error.message : 'Failed to generate molecule',
      })
    }
  },

  reset: () => {
    set(initialState)
  },

  setTool: (tool) => {
    set({ tool })
    // Auto-switch to add-atom when element selected
    if (tool === 'add-atom' && get().atomToAdd === null) {
      set({ atomToAdd: 'C' }) // Default to carbon
    }
  },

  resetTool: () => {
    set({ tool: 'select', atomToAdd: null })
  },

  setAtomToAdd: (element) => {
    set({ atomToAdd: element, tool: 'add-atom' })
  },

  setBondOrder: (order) => {
    set({ currentBondOrder: order })
  },

  setAutoBond: (value) => {
    set({ autoBond: value })
  },
}))

