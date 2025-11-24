import { describe, test, expect, beforeEach } from 'vitest'
import { useLabStore } from '../labStore'

describe('labStore', () => {
  beforeEach(() => {
    // Reset store before each test
    useLabStore.getState().resetMolecule()
  })

  test('add atom increases atom count and supports undo/redo', () => {
    const store = useLabStore.getState()
    store.resetMolecule()
    store.addAtom('C', [0, 0, 0])
    
    let mol = useLabStore.getState().molecule
    expect(mol.atoms.length).toBe(1)
    expect(mol.atoms[0].element).toBe('C')
    
    store.undo()
    mol = useLabStore.getState().molecule
    expect(mol.atoms.length).toBe(0)
    
    store.redo()
    mol = useLabStore.getState().molecule
    expect(mol.atoms.length).toBe(1)
  })

  test('addBond creates bond between atoms', () => {
    const store = useLabStore.getState()
    store.resetMolecule()
    const atom1 = store.addAtom('C', [0, 0, 0])
    const atom2 = store.addAtom('H', [1, 0, 0])
    
    const bond = store.addBond(atom1.id, atom2.id, 1)
    expect(bond).not.toBeNull()
    
    const mol = useLabStore.getState().molecule
    expect(mol.bonds.length).toBe(1)
    expect(mol.bonds[0].atom1).toBe(atom1.id)
    expect(mol.bonds[0].atom2).toBe(atom2.id)
  })

  test('deleteAtom removes atom and associated bonds', () => {
    const store = useLabStore.getState()
    store.resetMolecule()
    const atom1 = store.addAtom('C', [0, 0, 0])
    const atom2 = store.addAtom('H', [1, 0, 0])
    store.addBond(atom1.id, atom2.id, 1)
    
    store.deleteAtom(atom1.id)
    
    const mol = useLabStore.getState().molecule
    expect(mol.atoms.length).toBe(1)
    expect(mol.bonds.length).toBe(0)
  })
})

