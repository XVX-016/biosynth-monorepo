import { create } from 'zustand'
import { devtools } from 'zustand/middleware'
import { nanoid } from 'nanoid'
import type { Atom, Bond, Molecule, ToolName, Vec3 } from '../types/molecule'
import { canAddBond, validateValency } from '../chemcore/validation/valency'
import type { ReactionRule } from '../chemcore/reactions/simulator'
import { simulateReaction } from '../chemcore/reactions/simulator'

// Minimal helpers
const clone = <T,>(v: T) => JSON.parse(JSON.stringify(v)) as T

type LabState = {
  molecule: Molecule
  selectedAtomId: string | null
  selectedBondId: string | null
  currentTool: ToolName
  currentElement: string
  bondOrder: 1 | 2 | 3
  autoBond: boolean
  showGrid: boolean
  showValidation: boolean
  lastError: string | null

  // undo/redo
  undoStack: Molecule[]
  redoStack: Molecule[]

  // actions
  loadMolecule: (mol: Molecule) => void
  resetMolecule: () => void
  addAtom: (element: string, pos: Vec3) => Atom
  moveAtom: (atomId: string, pos: Vec3) => void
  deleteAtom: (atomId: string) => void
  addBond: (from: string, to: string, order?: 1 | 2 | 3) => Bond | null
  deleteBond: (bondId: string) => void
  applyReaction: (rule: ReactionRule) => void
  setTool: (t: ToolName) => void
  setBondOrder: (o: 1 | 2 | 3) => void
  setAutoBond: (v: boolean) => void
  setCurrentElement: (element: string) => void
  setSelectedAtomId: (id: string | null) => void
  setSelectedBondId: (id: string | null) => void
  setShowGrid: (v?: boolean) => void
  setShowValidation: (v?: boolean) => void
  clearError: () => void

  // undo/redo
  snapshot: () => void
  undo: () => void
  redo: () => void
}

const emptyMolecule = (): Molecule => ({ id: nanoid(), atoms: [], bonds: [], metadata: {} })

export const useLabStore = create<LabState>()(
  devtools((set, get) => ({
    molecule: emptyMolecule(),
    selectedAtomId: null,
    selectedBondId: null,
    currentTool: 'select',
    currentElement: 'C',
    bondOrder: 1,
    autoBond: true,
    showGrid: true,
    showValidation: true,
    lastError: null,
    undoStack: [],
    redoStack: [],

    loadMolecule: (mol) => {
      // eslint-disable-next-line no-console
      console.log('[LabStore] Loading molecule:', {
        id: mol.id,
        atomCount: mol.atoms.length,
        bondCount: mol.bonds.length,
        name: mol.metadata?.name,
        source: mol.metadata?.source,
      });
      // Validate molecule structure
      if (mol.atoms.length === 0 && mol.bonds.length > 0) {
        // eslint-disable-next-line no-console
        console.warn('[LabStore] Warning: molecule has bonds but no atoms');
      }
      // Validate atom positions
      mol.atoms.forEach((atom, idx) => {
        if (!atom.position || typeof atom.position.x !== 'number' || typeof atom.position.y !== 'number' || typeof atom.position.z !== 'number') {
          // eslint-disable-next-line no-console
          console.warn('[LabStore] Invalid atom position:', { atom, idx });
        }
      });
      set({ molecule: clone(mol), undoStack: [], redoStack: [] });
    },
    resetMolecule: () => set({ molecule: emptyMolecule(), undoStack: [], redoStack: [] }),

    snapshot: () => {
      const mol = get().molecule
      set((state) => ({ undoStack: [...state.undoStack, clone(mol)], redoStack: [] }))
    },

    undo: () => {
      const { undoStack, redoStack, molecule } = get()
      if (undoStack.length === 0) return
      const last = undoStack[undoStack.length - 1]
      set({ molecule: clone(last), undoStack: undoStack.slice(0, -1), redoStack: [...redoStack, clone(molecule)] })
    },

    redo: () => {
      const { undoStack, redoStack, molecule } = get()
      if (redoStack.length === 0) return
      const next = redoStack[redoStack.length - 1]
      set({ molecule: clone(next), redoStack: redoStack.slice(0, -1), undoStack: [...undoStack, clone(molecule)] })
    },

    clearError: () => set({ lastError: null }),

    addAtom: (element, pos) => {
      get().snapshot()
      const id = nanoid()
      const atom: Atom = { id, element, position: { x: pos[0], y: pos[1], z: pos[2] } }

      set((s) => ({ molecule: { ...s.molecule, atoms: [...s.molecule.atoms, atom] } }))

      // Legacy auto-bond logic removed
      return atom
    },

    moveAtom: (atomId, pos) => {
      get().snapshot()
      set((s) => ({
        molecule: {
          ...s.molecule,
          atoms: s.molecule.atoms.map(a => a.id === atomId ? { ...a, position: { x: pos[0], y: pos[1], z: pos[2] } } : a)
        }
      }))
    },

    deleteAtom: (atomId) => {
      get().snapshot()
      set((s) => ({
        molecule: {
          ...s.molecule,
          atoms: s.molecule.atoms.filter(a => a.id !== atomId),
          bonds: s.molecule.bonds.filter(b => b.from !== atomId && b.to !== atomId)
        }
      }))
    },

    addBond: (from, to, order) => {
      const state = get()
      const mol = state.molecule

      // 1. Sanity
      const exists = mol.bonds.find(x => (x.from === from && x.to === to) || (x.from === to && x.to === from))
      if (exists) return null

      // 2. Valency Validation
      const bondOrder = order || state.bondOrder || 1

      if (state.showValidation) {
        const validFrom = canAddBond(mol, from, bondOrder);
        const validTo = canAddBond(mol, to, bondOrder);

        if (!validFrom || !validTo) {
          set({ lastError: `Valency Limit Exceeded` });
          setTimeout(() => get().clearError(), 2000);
          return null;
        }
      }

      get().snapshot()
      const bond: Bond = { id: nanoid(), from, to, order: bondOrder }
      set((s) => ({ molecule: { ...s.molecule, bonds: [...s.molecule.bonds, bond] } }))
      return bond
    },

    deleteBond: (bondId) => {
      get().snapshot()
      set((s) => ({ molecule: { ...s.molecule, bonds: s.molecule.bonds.filter(b => b.id !== bondId) } }))
    },

    applyReaction: (rule) => {
      const state = get()
      const mol = state.molecule
      const result = simulateReaction(mol, rule)

      if (result) {
        const validation = validateValency(result)
        if (!validation.valid) {
          set({ lastError: `Reaction blocked: ${validation.reason || 'Valency limit'}` })
          setTimeout(() => get().clearError(), 3000)
          return
        }

        get().snapshot()
        set({ molecule: { ...mol, ...result } })
      } else {
        set({ lastError: "No applicable reaction sites found" })
        setTimeout(() => get().clearError(), 2000)
      }
    },

    setTool: (t) => set({ currentTool: t }),
    setBondOrder: (o) => set({ bondOrder: o }),
    setAutoBond: (v) => set({ autoBond: v }),
    setCurrentElement: (element) => set({ currentElement: element }),
    setSelectedAtomId: (id) => set({ selectedAtomId: id }),
    setSelectedBondId: (id) => set({ selectedBondId: id }),
    setShowGrid: (v) => set((s) => ({ showGrid: v ?? !s.showGrid })),
    setShowValidation: (v) => set((s) => ({ showValidation: v ?? !s.showValidation })),
  }), { name: 'LabStore' })
)
