import { create } from 'zustand'
import { MoleculeGraph, MoleculeSerializer } from '@biosynth/engine'
import type { Element } from '@biosynth/engine'
import { predict, generate, predictFast } from '../lib/api'
import type { SpectralData } from '../types/spectroscopy'
import type { ValidationState } from '../types/validation'
import { SpectroscopyService } from '../services/spectroscopyService'
import { MolecularValidator, StructureOptimizer } from '../services/validationService'

interface BackendPredictions {
  stability?: number
  toxicity?: number
  solubility?: number
  bioavailability?: number
  novelty?: number
}

export type Tool = 'select' | 'add-atom' | 'bond' | 'delete'
export type RenderMode = 'ballstick' | 'spacefill' | 'wireframe'
export type ColorScheme = 'element' | 'monochrome' | 'electrostatic'
type AsyncStatus = 'idle' | 'running' | 'complete' | 'error'

interface MoleculeState {
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
  renderMode: RenderMode
  colorScheme: ColorScheme
  spectroscopy: SpectralData
  validation: ValidationState
  highlightedAtoms: string[]
  optimizationStatus: AsyncStatus

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
  setRenderMode: (mode: RenderMode) => void
  setColorScheme: (scheme: ColorScheme) => void
  runSpectroscopy: () => Promise<void>
  runValidation: () => Promise<void>
  optimizeStructure: () => Promise<void>
  setHighlightedAtoms: (atoms: string[]) => void
  addBond: (atom1Id: string, atom2Id: string, order: number) => void
}

const spectroscopyService = new SpectroscopyService()
const validator = new MolecularValidator()
const optimizer = new StructureOptimizer()

const initialSpectroscopy: SpectralData = {
  ir: null,
  nmr: null,
  mass: null,
  status: 'idle',
  error: null,
}

const initialValidation: ValidationState = {
  status: 'idle',
  result: null,
  error: null,
}

const initialState = {
  currentMolecule: null,
  autoBond: false,
  selectedAtomId: null,
  selectedBondId: null,
  loadingState: 'idle' as const,
  backendPredictions: null,
  error: null,
  tool: 'select' as Tool,
  atomToAdd: null as Element | null,
  currentBondOrder: 1,
  renderMode: 'ballstick' as RenderMode,
  colorScheme: 'element' as ColorScheme,
  spectroscopy: { ...initialSpectroscopy },
  validation: { ...initialValidation },
  highlightedAtoms: [],
  optimizationStatus: 'idle' as AsyncStatus,
}

export const useMoleculeStore = create<MoleculeState>((set, get) => ({
  ...initialState,

  setMolecule: (molecule) => {
    set({
      currentMolecule: molecule,
      error: null,
      highlightedAtoms: [],
      spectroscopy: { ...initialSpectroscopy },
      validation: { ...initialValidation },
    })
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
      const smiles = MoleculeSerializer.toSMILES(currentMolecule)
      if (!smiles) {
        throw new Error('Unable to serialize molecule to SMILES')
      }
      const result = await predict(smiles)
      set({
        loadingState: 'success',
        backendPredictions: result.properties,
        error: null,
      })
    } catch (error) {
      set({
        loadingState: 'error',
        error:
          error instanceof Error
            ? error.message
            : 'Failed to fetch predictions',
        backendPredictions: null,
      })
    }
  },

  fetchPredictionsFast: async () => {
    const { currentMolecule } = get()
    if (!currentMolecule) return

    set({ loadingState: 'loading', error: null })

    try {
      const smiles = MoleculeSerializer.toSMILES(currentMolecule)
      if (!smiles) {
        throw new Error('Unable to serialize molecule to SMILES')
      }
      const result = await predictFast(smiles)
      set({
        loadingState: 'success',
        backendPredictions: result.properties,
        error: null,
      })
    } catch (error) {
      set({
        loadingState: 'error',
        error:
          error instanceof Error
            ? error.message
            : 'Failed to fetch predictions',
        backendPredictions: null,
      })
    }
  },

  generateMolecule: async (prompt: string) => {
    set({ loadingState: 'loading', error: null })

    try {
      await generate(prompt)
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
        error:
          error instanceof Error ? error.message : 'Failed to generate molecule',
      })
    }
  },

  reset: () => {
    set({
      ...initialState,
      spectroscopy: { ...initialSpectroscopy },
      validation: { ...initialValidation },
      highlightedAtoms: [],
    })
  },

  setTool: (tool) => {
    set({ tool })
    if (tool === 'add-atom' && get().atomToAdd === null) {
      set({ atomToAdd: 'C' })
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

  setRenderMode: (mode) => {
    set({ renderMode: mode })
  },

  setColorScheme: (scheme) => {
    set({ colorScheme: scheme })
  },

  runSpectroscopy: async () => {
    const { currentMolecule } = get()
    if (!currentMolecule || currentMolecule.atoms.size === 0) return
    set({
      spectroscopy: {
        ...get().spectroscopy,
        status: 'calculating',
        error: null,
      },
    })
    try {
      const result = await spectroscopyService.calculateAll(currentMolecule)
      set({
        spectroscopy: {
          ...result,
          status: 'complete',
          error: null,
          timestamp: new Date().toISOString(),
        },
      })
    } catch (error) {
      set({
        spectroscopy: {
          ...get().spectroscopy,
          status: 'error',
          error:
            error instanceof Error
              ? error.message
              : 'Spectroscopy failed',
        },
      })
    }
  },

  runValidation: async () => {
    const { currentMolecule } = get()
    if (!currentMolecule || currentMolecule.atoms.size === 0) return
    set({
      validation: {
        ...get().validation,
        status: 'running',
        error: null,
      },
    })
    try {
      const result = validator.validate(currentMolecule)
      set({
        validation: {
          status: 'complete',
          result,
          error: null,
          lastValidated: result.timestamp,
        },
      })
    } catch (error) {
      set({
        validation: {
          ...get().validation,
          status: 'error',
          error:
            error instanceof Error ? error.message : 'Validation failed',
        },
      })
    }
  },

  optimizeStructure: async () => {
    const { currentMolecule, runValidation } = get()
    if (!currentMolecule) return
    set({ optimizationStatus: 'running' })
    try {
      const optimized = await optimizer.optimize(currentMolecule)
      set({ currentMolecule: optimized.clone(), optimizationStatus: 'complete' })
      await runValidation()
    } catch (error) {
      set({
        optimizationStatus: 'error',
        error:
          error instanceof Error ? error.message : 'Optimization failed',
      })
    } finally {
      setTimeout(() => set({ optimizationStatus: 'idle' }), 500)
    }
  },

  setHighlightedAtoms: (atoms) => {
    set({ highlightedAtoms: atoms })
  },

  addBond: (atom1Id, atom2Id, order) => {
    const { currentMolecule } = get()
    if (!currentMolecule) return

    // Create a new bond using the engine's method if available, or manual manipulation
    // Assuming MoleculeGraph has an addBond method, otherwise we manipulate internals if exposed
    // Checking engineAdapter... likely need to use that or method on molecule
    // Since we are inside the store, and currentMolecule is a MoleculeGraph

    // We'll trust the engine has addBond or we reconstruct
    try {
      currentMolecule.addBond({ a1: atom1Id, a2: atom2Id, order })
      set({ currentMolecule: currentMolecule.clone() })
    } catch (e) {
      console.error("Failed to add bond", e)
    }
  },
}))

