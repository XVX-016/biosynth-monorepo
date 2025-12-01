/**
 * Molecule class - represents a complete molecule graph
 * 
 * Centralized molecule state management with validation and operations.
 */

import type { MoleculeState, Atom, Bond, ValidationError, ValidationResult } from './types'
import { AtomImpl } from './Atom'
import { BondImpl } from './Bond'
import { validateMolecule } from './validation/Validator'

export class Molecule {
  private atoms: Map<string, AtomImpl>
  private bonds: Map<string, BondImpl>
  private metadata: Record<string, any>

  constructor(state?: MoleculeState) {
    this.atoms = new Map()
    this.bonds = new Map()
    this.metadata = state?.metadata ?? {}

    if (state) {
      // Load atoms
      state.atoms.forEach((atomData, id) => {
        this.atoms.set(id, new AtomImpl(atomData))
      })

      // Load bonds
      state.bonds.forEach((bondData, id) => {
        this.bonds.set(id, new BondImpl(bondData))
      })
    }
  }

  /**
   * Get all atoms
   */
  getAtoms(): AtomImpl[] {
    return Array.from(this.atoms.values())
  }

  /**
   * Get atom by ID
   */
  getAtom(id: string): AtomImpl | null {
    return this.atoms.get(id) ?? null
  }

  /**
   * Get all bonds
   */
  getBonds(): BondImpl[] {
    return Array.from(this.bonds.values())
  }

  /**
   * Get bond by ID
   */
  getBond(id: string): BondImpl | null {
    return this.bonds.get(id) ?? null
  }

  /**
   * Get bonds connected to an atom
   */
  getBondsForAtom(atomId: string): BondImpl[] {
    return Array.from(this.bonds.values()).filter(bond => bond.connects(atomId))
  }

  /**
   * Get neighbors of an atom
   */
  getNeighbors(atomId: string): AtomImpl[] {
    const neighbors: AtomImpl[] = []
    this.getBondsForAtom(atomId).forEach(bond => {
      const otherId = bond.getOtherAtom(atomId)
      if (otherId) {
        const atom = this.getAtom(otherId)
        if (atom) neighbors.push(atom)
      }
    })
    return neighbors
  }

  /**
   * Add an atom
   */
  addAtom(atom: Atom): AtomImpl {
    const atomImpl = new AtomImpl(atom)
    this.atoms.set(atom.id, atomImpl)
    return atomImpl
  }

  /**
   * Remove an atom and all its bonds
   */
  removeAtom(atomId: string): void {
    // Remove all bonds connected to this atom
    const bondsToRemove = this.getBondsForAtom(atomId)
    bondsToRemove.forEach(bond => {
      this.bonds.delete(bond.id)
    })

    // Remove the atom
    this.atoms.delete(atomId)
  }

  /**
   * Add a bond
   */
  addBond(bond: Bond): BondImpl {
    // Validate atoms exist
    if (!this.atoms.has(bond.atom1)) {
      throw new Error(`Atom ${bond.atom1} does not exist`)
    }
    if (!this.atoms.has(bond.atom2)) {
      throw new Error(`Atom ${bond.atom2} does not exist`)
    }

    const bondImpl = new BondImpl(bond)
    this.bonds.set(bond.id, bondImpl)
    return bondImpl
  }

  /**
   * Remove a bond
   */
  removeBond(bondId: string): void {
    this.bonds.delete(bondId)
  }

  /**
   * Update atom position
   */
  updateAtomPosition(atomId: string, position: [number, number, number]): void {
    const atom = this.atoms.get(atomId)
    if (!atom) {
      throw new Error(`Atom ${atomId} does not exist`)
    }
    const updated = atom.update({ position })
    this.atoms.set(atomId, updated)
  }

  /**
   * Update atom properties (immutable)
   */
  updateAtom(atomId: string, updates: Partial<Atom>): Molecule {
    const atom = this.atoms.get(atomId)
    if (!atom) return this

    const updatedAtom = new AtomImpl({
      ...atom.toJSON(),
      ...updates,
    })

    const newMolecule = new Molecule()
    this.atoms.forEach((a, id) => {
      newMolecule.atoms.set(id, id === atomId ? updatedAtom : a)
    })
    this.bonds.forEach((bond, id) => {
      newMolecule.bonds.set(id, bond)
    })
    newMolecule.metadata = { ...this.metadata }

    return newMolecule
  }

  /**
   * Update bond properties (immutable)
   */
  updateBond(bondId: string, updates: Partial<Bond>): Molecule {
    const bond = this.bonds.get(bondId)
    if (!bond) return this

    const updatedBond = new BondImpl({
      ...bond.toJSON(),
      ...updates,
    })

    const newMolecule = new Molecule()
    this.atoms.forEach((atom, id) => {
      newMolecule.atoms.set(id, atom)
    })
    this.bonds.forEach((b, id) => {
      newMolecule.bonds.set(id, id === bondId ? updatedBond : b)
    })
    newMolecule.metadata = { ...this.metadata }

    return newMolecule
  }

  /**
   * Clear all atoms and bonds
   */
  clear(): void {
    this.atoms.clear()
    this.bonds.clear()
    this.metadata = {}
  }

  /**
   * Get molecule size
   */
  get size(): { atoms: number; bonds: number } {
    return {
      atoms: this.atoms.size,
      bonds: this.bonds.size,
    }
  }

  /**
   * Check if molecule is empty
   */
  isEmpty(): boolean {
    return this.atoms.size === 0
  }

  /**
   * Validate molecule structure
   */
  validate(): ValidationResult {
    return validateMolecule(this)
  }

  /**
   * Get metadata
   */
  getMetadata(): Record<string, any> {
    return { ...this.metadata }
  }

  /**
   * Set metadata
   */
  setMetadata(metadata: Record<string, any>): void {
    this.metadata = { ...this.metadata, ...metadata }
  }

  /**
   * Convert to plain state object
   */
  toState(): MoleculeState {
    const atoms = new Map<string, Atom>()
    const bonds = new Map<string, Bond>()

    this.atoms.forEach((atom, id) => {
      atoms.set(id, atom.toJSON())
    })

    this.bonds.forEach((bond, id) => {
      bonds.set(id, bond.toJSON())
    })

    return {
      atoms,
      bonds,
      metadata: { ...this.metadata },
    }
  }

  /**
   * Create a copy of the molecule
   */
  clone(): Molecule {
    return new Molecule(this.toState())
  }

  /**
   * Generate 2D layout (async, uses backend)
   * 
   * This is a convenience method that calls the layout service.
   * For better control, use generate2DLayout from layout service directly.
   */
  async autoLayout(options?: { method?: 'coordgen' | 'rdkit'; spacing?: number }): Promise<Molecule> {
    // Dynamic import to avoid circular dependency
    const { generate2DLayout } = await import('./layout')
    return generate2DLayout(this, options)
  }
}

