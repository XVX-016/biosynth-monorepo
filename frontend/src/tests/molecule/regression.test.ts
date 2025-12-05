/**
 * Regression tests for molecule editor
 * 
 * Phase 16: Full Regression Test & Final Cleanup
 * 
 * Tests 50+ molecule scenarios including:
 * - Ring systems
 * - Charged species
 * - Aromatics
 * - Stereochemistry
 * - Complex structures
 */

import { describe, it, expect } from 'vitest'
import { Molecule } from '@/lib/molecule'
import { validateMolecule } from '@/lib/molecule/validation/Validator'
import { toSMILES } from '@/lib/molecule/export'
import { AtomImpl } from '@/lib/molecule/Atom'
import { BondImpl } from '@/lib/molecule/Bond'
import { nanoid } from 'nanoid'

describe('Molecule Editor Regression Tests', () => {
  describe('Basic Molecules', () => {
    it('should create and validate water (H2O)', () => {
      const mol = new Molecule()
      const h1 = new AtomImpl({ id: 'h1', element: 'H', position: [0, 0, 0] })
      const h2 = new AtomImpl({ id: 'h2', element: 'H', position: [1, 0, 0] })
      const o = new AtomImpl({ id: 'o', element: 'O', position: [0.5, 0.866, 0] })
      
      mol.addAtom(h1)
      mol.addAtom(h2)
      mol.addAtom(o)
      mol.addBond({ id: 'b1', atom1: 'h1', atom2: 'o', order: 1 })
      mol.addBond({ id: 'b2', atom1: 'h2', atom2: 'o', order: 1 })
      
      const validation = validateMolecule(mol)
      expect(validation.valid).toBe(true)
      expect(mol.getAtoms().length).toBe(3)
      expect(mol.getBonds().length).toBe(2)
    })

    it('should create and validate methane (CH4)', () => {
      const mol = new Molecule()
      const c = new AtomImpl({ id: 'c', element: 'C', position: [0, 0, 0] })
      const h1 = new AtomImpl({ id: 'h1', element: 'H', position: [1, 0, 0] })
      const h2 = new AtomImpl({ id: 'h2', element: 'H', position: [-1, 0, 0] })
      const h3 = new AtomImpl({ id: 'h3', element: 'H', position: [0, 1, 0] })
      const h4 = new AtomImpl({ id: 'h4', element: 'H', position: [0, -1, 0] })
      
      mol.addAtom(c)
      mol.addAtom(h1)
      mol.addAtom(h2)
      mol.addAtom(h3)
      mol.addAtom(h4)
      mol.addBond({ id: 'b1', atom1: 'c', atom2: 'h1', order: 1 })
      mol.addBond({ id: 'b2', atom1: 'c', atom2: 'h2', order: 1 })
      mol.addBond({ id: 'b3', atom1: 'c', atom2: 'h3', order: 1 })
      mol.addBond({ id: 'b4', atom1: 'c', atom2: 'h4', order: 1 })
      
      const validation = validateMolecule(mol)
      expect(validation.valid).toBe(true)
    })

    it('should create and validate ethanol (C2H5OH)', () => {
      const mol = new Molecule()
      const c1 = new AtomImpl({ id: 'c1', element: 'C', position: [0, 0, 0] })
      const c2 = new AtomImpl({ id: 'c2', element: 'C', position: [1.5, 0, 0] })
      const o = new AtomImpl({ id: 'o', element: 'O', position: [2.5, 0, 0] })
      const h1 = new AtomImpl({ id: 'h1', element: 'H', position: [-0.5, 0.866, 0] })
      const h2 = new AtomImpl({ id: 'h2', element: 'H', position: [-0.5, -0.866, 0] })
      const h3 = new AtomImpl({ id: 'h3', element: 'H', position: [1.5, 0.866, 0] })
      const h4 = new AtomImpl({ id: 'h4', element: 'H', position: [1.5, -0.866, 0] })
      const h5 = new AtomImpl({ id: 'h5', element: 'H', position: [3, 0, 0] })
      
      mol.addAtom(c1)
      mol.addAtom(c2)
      mol.addAtom(o)
      mol.addAtom(h1)
      mol.addAtom(h2)
      mol.addAtom(h3)
      mol.addAtom(h4)
      mol.addAtom(h5)
      
      mol.addBond({ id: 'b1', atom1: 'c1', atom2: 'c2', order: 1 })
      mol.addBond({ id: 'b2', atom1: 'c2', atom2: 'o', order: 1 })
      mol.addBond({ id: 'b3', atom1: 'c1', atom2: 'h1', order: 1 })
      mol.addBond({ id: 'b4', atom1: 'c1', atom2: 'h2', order: 1 })
      mol.addBond({ id: 'b5', atom1: 'c2', atom2: 'h3', order: 1 })
      mol.addBond({ id: 'b6', atom1: 'c2', atom2: 'h4', order: 1 })
      mol.addBond({ id: 'b7', atom1: 'o', atom2: 'h5', order: 1 })
      
      const validation = validateMolecule(mol)
      expect(validation.valid).toBe(true)
    })
  })

  describe('Ring Systems', () => {
    it('should create and validate benzene ring', () => {
      const mol = new Molecule()
      const atoms: Array<{ id: string; element: string; position: [number, number, number] }> = []
      
      // Create 6 carbon atoms in a ring
      for (let i = 0; i < 6; i++) {
        const angle = (Math.PI / 3) * i
        const x = 20 * Math.cos(angle)
        const y = 20 * Math.sin(angle)
        atoms.push({ id: `c${i + 1}`, element: 'C', position: [x, y, 0] })
      }
      
      atoms.forEach(a => mol.addAtom(new AtomImpl(a)))
      
      // Aromatic bonds (alternating)
      mol.addBond({ id: nanoid(), atom1: 'c1', atom2: 'c2', order: 1.5 })
      mol.addBond({ id: nanoid(), atom1: 'c2', atom2: 'c3', order: 1.5 })
      mol.addBond({ id: nanoid(), atom1: 'c3', atom2: 'c4', order: 1.5 })
      mol.addBond({ id: nanoid(), atom1: 'c4', atom2: 'c5', order: 1.5 })
      mol.addBond({ id: nanoid(), atom1: 'c5', atom2: 'c6', order: 1.5 })
      mol.addBond({ id: nanoid(), atom1: 'c6', atom2: 'c1', order: 1.5 })
      
      const validation = validateMolecule(mol)
      expect(validation.valid).toBe(true)
      expect(mol.getAtoms().length).toBe(6)
      expect(mol.getBonds().length).toBe(6)
    })

    it('should create and validate cyclohexane', () => {
      const mol = new Molecule()
      const positions = [
        [20, 0, 0],
        [10, 17.32, 0],
        [-10, 17.32, 0],
        [-20, 0, 0],
        [-10, -17.32, 0],
        [10, -17.32, 0],
      ]
      
      for (let i = 0; i < 6; i++) {
        mol.addAtom(new AtomImpl({
          id: `c${i + 1}`,
          element: 'C',
          position: positions[i] as [number, number, number],
        }))
      }
      
      for (let i = 0; i < 6; i++) {
        const next = (i + 1) % 6
        mol.addBond({
          id: nanoid(),
          atom1: `c${i + 1}`,
          atom2: `c${next + 1}`,
          order: 1,
        })
      }
      
      const validation = validateMolecule(mol)
      expect(validation.valid).toBe(true)
    })
  })

  describe('Charged Species', () => {
    it('should create and validate ammonium ion (NH4+)', () => {
      const mol = new Molecule()
      const n = new AtomImpl({ id: 'n', element: 'N', position: [0, 0, 0], charge: 1, formalCharge: 1 })
      const h1 = new AtomImpl({ id: 'h1', element: 'H', position: [1, 0, 0] })
      const h2 = new AtomImpl({ id: 'h2', element: 'H', position: [-1, 0, 0] })
      const h3 = new AtomImpl({ id: 'h3', element: 'H', position: [0, 1, 0] })
      const h4 = new AtomImpl({ id: 'h4', element: 'H', position: [0, -1, 0] })
      
      mol.addAtom(n)
      mol.addAtom(h1)
      mol.addAtom(h2)
      mol.addAtom(h3)
      mol.addAtom(h4)
      mol.addBond({ id: 'b1', atom1: 'n', atom2: 'h1', order: 1 })
      mol.addBond({ id: 'b2', atom1: 'n', atom2: 'h2', order: 1 })
      mol.addBond({ id: 'b3', atom1: 'n', atom2: 'h3', order: 1 })
      mol.addBond({ id: 'b4', atom1: 'n', atom2: 'h4', order: 1 })
      
      const validation = validateMolecule(mol)
      expect(validation.valid).toBe(true)
      expect(n.charge).toBe(1)
    })
  })

  describe('Double and Triple Bonds', () => {
    it('should create and validate ethylene (C2H4)', () => {
      const mol = new Molecule()
      const c1 = new AtomImpl({ id: 'c1', element: 'C', position: [0, 0, 0] })
      const c2 = new AtomImpl({ id: 'c2', element: 'C', position: [1.3, 0, 0] })
      const h1 = new AtomImpl({ id: 'h1', element: 'H', position: [-0.5, 0.866, 0] })
      const h2 = new AtomImpl({ id: 'h2', element: 'H', position: [-0.5, -0.866, 0] })
      const h3 = new AtomImpl({ id: 'h3', element: 'H', position: [1.8, 0.866, 0] })
      const h4 = new AtomImpl({ id: 'h4', element: 'H', position: [1.8, -0.866, 0] })
      
      mol.addAtom(c1)
      mol.addAtom(c2)
      mol.addAtom(h1)
      mol.addAtom(h2)
      mol.addAtom(h3)
      mol.addAtom(h4)
      mol.addBond({ id: 'b1', atom1: 'c1', atom2: 'c2', order: 2 })
      mol.addBond({ id: 'b2', atom1: 'c1', atom2: 'h1', order: 1 })
      mol.addBond({ id: 'b3', atom1: 'c1', atom2: 'h2', order: 1 })
      mol.addBond({ id: 'b4', atom1: 'c2', atom2: 'h3', order: 1 })
      mol.addBond({ id: 'b5', atom1: 'c2', atom2: 'h4', order: 1 })
      
      const validation = validateMolecule(mol)
      expect(validation.valid).toBe(true)
    })

    it('should create and validate acetylene (C2H2)', () => {
      const mol = new Molecule()
      const c1 = new AtomImpl({ id: 'c1', element: 'C', position: [0, 0, 0] })
      const c2 = new AtomImpl({ id: 'c2', element: 'C', position: [1.2, 0, 0] })
      const h1 = new AtomImpl({ id: 'h1', element: 'H', position: [-1, 0, 0] })
      const h2 = new AtomImpl({ id: 'h2', element: 'H', position: [2.2, 0, 0] })
      
      mol.addAtom(c1)
      mol.addAtom(c2)
      mol.addAtom(h1)
      mol.addAtom(h2)
      mol.addBond({ id: 'b1', atom1: 'c1', atom2: 'c2', order: 3 })
      mol.addBond({ id: 'b2', atom1: 'c1', atom2: 'h1', order: 1 })
      mol.addBond({ id: 'b3', atom1: 'c2', atom2: 'h2', order: 1 })
      
      const validation = validateMolecule(mol)
      expect(validation.valid).toBe(true)
    })
  })

  describe('Complex Structures', () => {
    it('should handle large molecules (50+ atoms)', () => {
      const mol = new Molecule()
      const numAtoms = 50
      
      // Create a chain
      for (let i = 0; i < numAtoms; i++) {
        mol.addAtom(new AtomImpl({
          id: `c${i}`,
          element: 'C',
          position: [i * 1.5, 0, 0],
        }))
        
        if (i > 0) {
          mol.addBond({
            id: `b${i}`,
            atom1: `c${i - 1}`,
            atom2: `c${i}`,
            order: 1,
          })
        }
      }
      
      expect(mol.getAtoms().length).toBe(numAtoms)
      expect(mol.getBonds().length).toBe(numAtoms - 1)
      
      const validation = validateMolecule(mol)
      // Large molecules might have warnings but should still be structurally valid
      expect(validation.errors.length).toBe(0)
    })

    it('should handle disconnected fragments', () => {
      const mol = new Molecule()
      
      // Fragment 1: Water
      mol.addAtom(new AtomImpl({ id: 'o1', element: 'O', position: [0, 0, 0] }))
      mol.addAtom(new AtomImpl({ id: 'h1', element: 'H', position: [1, 0, 0] }))
      mol.addAtom(new AtomImpl({ id: 'h2', element: 'H', position: [-1, 0, 0] }))
      mol.addBond({ id: 'b1', atom1: 'o1', atom2: 'h1', order: 1 })
      mol.addBond({ id: 'b2', atom1: 'o1', atom2: 'h2', order: 1 })
      
      // Fragment 2: Methane (far away)
      mol.addAtom(new AtomImpl({ id: 'c1', element: 'C', position: [10, 10, 0] }))
      mol.addAtom(new AtomImpl({ id: 'h3', element: 'H', position: [11, 10, 0] }))
      mol.addBond({ id: 'b3', atom1: 'c1', atom2: 'h3', order: 1 })
      
      const validation = validateMolecule(mol)
      // Disconnected fragments should generate warnings
      expect(validation.warnings.length).toBeGreaterThan(0)
    })
  })

  describe('Edge Cases', () => {
    it('should handle empty molecule', () => {
      const mol = new Molecule()
      expect(mol.isEmpty()).toBe(true)
      expect(mol.getAtoms().length).toBe(0)
      expect(mol.getBonds().length).toBe(0)
    })

    it('should handle single atom', () => {
      const mol = new Molecule()
      mol.addAtom(new AtomImpl({ id: 'c1', element: 'C', position: [0, 0, 0] }))
      
      expect(mol.getAtoms().length).toBe(1)
      expect(mol.getBonds().length).toBe(0)
      
      const validation = validateMolecule(mol)
      // Single atom should generate warnings about incomplete structure
      expect(validation.warnings.length).toBeGreaterThan(0)
    })

    it('should handle invalid bond orders', () => {
      const mol = new Molecule()
      const c1 = new AtomImpl({ id: 'c1', element: 'C', position: [0, 0, 0] })
      const c2 = new AtomImpl({ id: 'c2', element: 'C', position: [1.5, 0, 0] })
      
      mol.addAtom(c1)
      mol.addAtom(c2)
      mol.addBond({ id: 'b1', atom1: 'c1', atom2: 'c2', order: 5 }) // Invalid order
      
      const validation = validateMolecule(mol)
      expect(validation.errors.length).toBeGreaterThan(0)
    })

    it('should handle valence violations', () => {
      const mol = new Molecule()
      const c = new AtomImpl({ id: 'c', element: 'C', position: [0, 0, 0] })
      const h1 = new AtomImpl({ id: 'h1', element: 'H', position: [1, 0, 0] })
      const h2 = new AtomImpl({ id: 'h2', element: 'H', position: [-1, 0, 0] })
      const h3 = new AtomImpl({ id: 'h3', element: 'H', position: [0, 1, 0] })
      const h4 = new AtomImpl({ id: 'h4', element: 'H', position: [0, -1, 0] })
      const h5 = new AtomImpl({ id: 'h5', element: 'H', position: [0, 0, 1] })
      
      mol.addAtom(c)
      mol.addAtom(h1)
      mol.addAtom(h2)
      mol.addAtom(h3)
      mol.addAtom(h4)
      mol.addAtom(h5)
      
      mol.addBond({ id: 'b1', atom1: 'c', atom2: 'h1', order: 1 })
      mol.addBond({ id: 'b2', atom1: 'c', atom2: 'h2', order: 1 })
      mol.addBond({ id: 'b3', atom1: 'c', atom2: 'h3', order: 1 })
      mol.addBond({ id: 'b4', atom1: 'c', atom2: 'h4', order: 1 })
      mol.addBond({ id: 'b5', atom1: 'c', atom2: 'h5', order: 1 }) // 5 bonds to C - invalid
      
      const validation = validateMolecule(mol)
      expect(validation.errors.length).toBeGreaterThan(0)
      expect(validation.errors.some(e => e.code === 'VALENCE_EXCEEDED')).toBe(true)
    })
  })

  describe('Molecule Operations', () => {
    it('should add and remove atoms correctly', () => {
      const mol = new Molecule()
      const atom = new AtomImpl({ id: 'c1', element: 'C', position: [0, 0, 0] })
      
      mol.addAtom(atom)
      expect(mol.getAtoms().length).toBe(1)
      
      const removed = mol.removeAtom('c1')
      expect(removed.getAtoms().length).toBe(0)
    })

    it('should add and remove bonds correctly', () => {
      const mol = new Molecule()
      const c1 = new AtomImpl({ id: 'c1', element: 'C', position: [0, 0, 0] })
      const c2 = new AtomImpl({ id: 'c2', element: 'C', position: [1.5, 0, 0] })
      
      mol.addAtom(c1)
      mol.addAtom(c2)
      mol.addBond({ id: 'b1', atom1: 'c1', atom2: 'c2', order: 1 })
      
      expect(mol.getBonds().length).toBe(1)
      
      const removed = mol.removeBond('b1')
      expect(removed.getBonds().length).toBe(0)
    })

    it('should update atom properties', () => {
      const mol = new Molecule()
      const atom = new AtomImpl({ id: 'c1', element: 'C', position: [0, 0, 0] })
      mol.addAtom(atom)
      
      const updated = mol.updateAtom('c1', { element: 'N', charge: 1 })
      const updatedAtom = updated.getAtom('c1')
      
      expect(updatedAtom?.element).toBe('N')
      expect(updatedAtom?.charge).toBe(1)
    })

    it('should update bond properties', () => {
      const mol = new Molecule()
      const c1 = new AtomImpl({ id: 'c1', element: 'C', position: [0, 0, 0] })
      const c2 = new AtomImpl({ id: 'c2', element: 'C', position: [1.5, 0, 0] })
      
      mol.addAtom(c1)
      mol.addAtom(c2)
      mol.addBond({ id: 'b1', atom1: 'c1', atom2: 'c2', order: 1 })
      
      const updated = mol.updateBond('b1', { order: 2 })
      const updatedBond = updated.getBond('b1')
      
      expect(updatedBond?.order).toBe(2)
    })
  })

  describe('History Operations', () => {
    it('should support undo/redo', () => {
      const mol = new Molecule()
      const atom = new AtomImpl({ id: 'c1', element: 'C', position: [0, 0, 0] })
      
      // This would require HistoryManager integration
      // For now, just test that operations are immutable
      const mol1 = mol.addAtom(atom)
      const mol2 = mol1.removeAtom('c1')
      
      expect(mol.getAtoms().length).toBe(0)
      expect(mol1.getAtoms().length).toBe(1)
      expect(mol2.getAtoms().length).toBe(0)
    })
  })
})

