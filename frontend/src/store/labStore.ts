import { create } from 'zustand'
import { devtools } from 'zustand/middleware'
import { nanoid } from 'nanoid'
import type { Atom, Bond, Molecule, ToolName, Vec3 } from '../types/molecule'
import { computeAutoBonds } from '../utils/bondingEngine'

// Minimal helpers
const clone = <T,>(v: T) => JSON.parse(JSON.stringify(v)) as T

type LabState = {
  molecule: Molecule
  selectedAtomId: string | null
  selectedBondId: string | null
  currentTool: ToolName
  currentElement: string
  autoBond: boolean
  // undo/redo
  undoStack: Molecule[]
  redoStack: Molecule[]

  // actions
  loadMolecule: (mol: Molecule) => void
  resetMolecule: () => void
  addAtom: (element: string, pos: Vec3) => Atom
  moveAtom: (atomId: string, pos: Vec3) => void
  deleteAtom: (atomId: string) => void
  addBond: (a: string, b: string, order?: 1|2|3) => Bond | null
  deleteBond: (bondId: string) => void
  setTool: (t: ToolName) => void
  setAutoBond: (v: boolean) => void
  setCurrentElement: (element: string) => void
  setSelectedAtomId: (id: string | null) => void
  setSelectedBondId: (id: string | null) => void

  // undo/redo
  snapshot: () => void
  undo: () => void
  redo: () => void
}

const emptyMolecule = (): Molecule => ({ atoms: [], bonds: [], metadata: {} })

export const useLabStore = create<LabState>()(
  devtools((set, get) => ({
    molecule: emptyMolecule(),
    selectedAtomId: null,
    selectedBondId: null,
    currentTool: 'select',
    currentElement: 'C',
    autoBond: true,
    undoStack: [],
    redoStack: [],

    loadMolecule: (mol) => set({ molecule: clone(mol), undoStack: [], redoStack: [] }),
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

    addAtom: (element, pos) => {
      get().snapshot()
      const id = nanoid()
      const atom: Atom = { id, element, position: pos, charge: 0 }
      set((s) => ({ molecule: { ...s.molecule, atoms: [...s.molecule.atoms, atom] } }))
      
      // Auto-bond if enabled
      const state = get()
      if (state.autoBond) {
        const mol = state.molecule
        const candidates = computeAutoBonds(mol)
        candidates.forEach(c => {
          const exists = mol.bonds.find(b => 
            (b.atom1 === c.a && b.atom2 === c.b) || (b.atom1 === c.b && b.atom2 === c.a)
          )
          if (!exists) {
            const bid = nanoid()
            set((s) => ({ 
              molecule: { 
                ...s.molecule, 
                bonds: [...s.molecule.bonds, { id: bid, atom1: c.a, atom2: c.b, order: 1 }] 
              } 
            }))
          }
        })
      }
      
      return atom
    },

    moveAtom: (atomId, pos) => {
      get().snapshot()
      set((s) => ({
        molecule: {
          ...s.molecule,
          atoms: s.molecule.atoms.map(a => a.id === atomId ? { ...a, position: pos } : a)
        }
      }))
    },

    deleteAtom: (atomId) => {
      get().snapshot()
      set((s) => ({
        molecule: {
          ...s.molecule,
          atoms: s.molecule.atoms.filter(a => a.id !== atomId),
          bonds: s.molecule.bonds.filter(b => b.atom1 !== atomId && b.atom2 !== atomId)
        }
      }))
    },

    addBond: (a, b, order = 1) => {
      // simple sanity: no dup bonds
      const mol = get().molecule
      const exists = mol.bonds.find(x => (x.atom1 === a && x.atom2 === b) || (x.atom1 === b && x.atom2 === a))
      if (exists) return null
      get().snapshot()
      const bond: Bond = { id: nanoid(), atom1: a, atom2: b, order }
      set((s) => ({ molecule: { ...s.molecule, bonds: [...s.molecule.bonds, bond] } }))
      return bond
    },

    deleteBond: (bondId) => {
      get().snapshot()
      set((s) => ({ molecule: { ...s.molecule, bonds: s.molecule.bonds.filter(b => b.id !== bondId) } }))
    },

    setTool: (t) => set({ currentTool: t }),
    setAutoBond: (v) => set({ autoBond: v }),
    setCurrentElement: (element) => set({ currentElement: element }),
    setSelectedAtomId: (id) => set({ selectedAtomId: id }),
    setSelectedBondId: (id) => set({ selectedBondId: id }),
  }), { name: 'LabStore' })
)

