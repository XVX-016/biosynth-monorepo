/**
 * Command - Base command interface for undo/redo
 * 
 * Phase 4: Undo/Redo System
 * 
 * All molecule operations are commands that can be executed and undone.
 */

import { Molecule } from '../Molecule'

export interface Command {
  /**
   * Execute the command
   */
  execute(molecule: Molecule): void

  /**
   * Undo the command
   */
  undo(molecule: Molecule): void

  /**
   * Get command description (for UI)
   */
  getDescription(): string
}

/**
 * Command to add an atom
 */
export class AddAtomCommand implements Command {
  constructor(
    private atomId: string,
    private element: string,
    private position: [number, number, number],
  ) {}

  execute(molecule: Molecule): void {
    molecule.addAtom({
      id: this.atomId,
      element: this.element,
      position: this.position,
    })
  }

  undo(molecule: Molecule): void {
    molecule.removeAtom(this.atomId)
  }

  getDescription(): string {
    return `Add ${this.element} atom`
  }
}

/**
 * Command to remove an atom
 */
export class RemoveAtomCommand implements Command {
  private atomData: any
  private bondsData: any[] = []

  constructor(private atomId: string) {}

  execute(molecule: Molecule): void {
    // Store atom data for undo
    const atom = molecule.getAtom(this.atomId)
    if (atom) {
      this.atomData = atom.toJSON()
    }

    // Store all bonds connected to this atom
    const bonds = molecule.getBondsForAtom(this.atomId)
    this.bondsData = bonds.map(bond => bond.toJSON())

    // Remove atom (this also removes all connected bonds)
    molecule.removeAtom(this.atomId)
  }

  undo(molecule: Molecule): void {
    // Restore atom
    if (this.atomData) {
      molecule.addAtom(this.atomData)
    }

    // Restore bonds
    this.bondsData.forEach(bondData => {
      molecule.addBond(bondData)
    })
  }

  getDescription(): string {
    return 'Remove atom'
  }
}

/**
 * Command to add a bond
 */
export class AddBondCommand implements Command {
  constructor(
    private bondId: string,
    private atom1Id: string,
    private atom2Id: string,
    private order: number,
  ) {}

  execute(molecule: Molecule): void {
    molecule.addBond({
      id: this.bondId,
      atom1: this.atom1Id,
      atom2: this.atom2Id,
      order: this.order,
    })
  }

  undo(molecule: Molecule): void {
    molecule.removeBond(this.bondId)
  }

  getDescription(): string {
    return `Add bond (order ${this.order})`
  }
}

/**
 * Command to remove a bond
 */
export class RemoveBondCommand implements Command {
  private bondData: any

  constructor(private bondId: string) {}

  execute(molecule: Molecule): void {
    // Store bond data for undo
    const bond = molecule.getBond(this.bondId)
    if (bond) {
      this.bondData = bond.toJSON()
    }

    molecule.removeBond(this.bondId)
  }

  undo(molecule: Molecule): void {
    if (this.bondData) {
      molecule.addBond(this.bondData)
    }
  }

  getDescription(): string {
    return 'Remove bond'
  }
}

/**
 * Command to move an atom
 */
export class MoveAtomCommand implements Command {
  private oldPosition: [number, number, number]

  constructor(
    private atomId: string,
    private newPosition: [number, number, number],
  ) {}

  execute(molecule: Molecule): void {
    const atom = molecule.getAtom(this.atomId)
    if (atom) {
      this.oldPosition = [...atom.position] as [number, number, number]
      molecule.updateAtomPosition(this.atomId, this.newPosition)
    }
  }

  undo(molecule: Molecule): void {
    if (this.oldPosition) {
      molecule.updateAtomPosition(this.atomId, this.oldPosition)
    }
  }

  getDescription(): string {
    return 'Move atom'
  }
}

/**
 * Command to update atom properties
 */
export class UpdateAtomCommand implements Command {
  private oldData: any

  constructor(
    private atomId: string,
    private updates: any,
  ) {}

  execute(molecule: Molecule): void {
    const atom = molecule.getAtom(this.atomId)
    if (atom) {
      this.oldData = atom.toJSON()
      molecule.updateAtom(this.atomId, this.updates)
    }
  }

  undo(molecule: Molecule): void {
    if (this.oldData) {
      molecule.updateAtom(this.atomId, this.oldData)
    }
  }

  getDescription(): string {
    return 'Update atom'
  }
}

/**
 * Command to update bond properties
 */
export class UpdateBondCommand implements Command {
  private oldData: any

  constructor(
    private bondId: string,
    private updates: any,
  ) {}

  execute(molecule: Molecule): void {
    const bond = molecule.getBond(this.bondId)
    if (bond) {
      this.oldData = bond.toJSON()
      molecule.updateBond(this.bondId, this.updates)
    }
  }

  undo(molecule: Molecule): void {
    if (this.oldData) {
      molecule.updateBond(this.bondId, this.oldData)
    }
  }

  getDescription(): string {
    return 'Update bond'
  }
}

/**
 * Command to clear molecule
 */
export class ClearMoleculeCommand implements Command {
  private savedState: any

  execute(molecule: Molecule): void {
    // Save current state
    this.savedState = molecule.toState()
    molecule.clear()
  }

  undo(molecule: Molecule): void {
    if (this.savedState) {
      // Restore from saved state
      const restored = new Molecule(this.savedState)
      molecule.clear()
      // Copy all atoms and bonds
      restored.getAtoms().forEach(atom => {
        molecule.addAtom(atom.toJSON())
      })
      restored.getBonds().forEach(bond => {
        molecule.addBond(bond.toJSON())
      })
    }
  }

  getDescription(): string {
    return 'Clear molecule'
  }
}

