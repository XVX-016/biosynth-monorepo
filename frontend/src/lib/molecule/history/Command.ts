/**
 * Command Pattern Implementation
 * 
 * Phase 4: Undo/Redo System
 * Phase 12: Added UpdateAtomCommand and UpdateBondCommand
 * 
 * Each command encapsulates an operation and its inverse.
 */

import type { Molecule } from '../Molecule'
import type { Atom, Bond } from '../types'
import { nanoid } from 'nanoid'

export interface Command {
  execute(molecule: Molecule): Molecule
  undo(molecule: Molecule): Molecule
  description: string
}

export class AddAtomCommand implements Command {
  private atom: Atom
  constructor(atom: Atom) {
    this.atom = atom
  }
  execute(molecule: Molecule): Molecule {
    return molecule.addAtom(this.atom)
  }
  undo(molecule: Molecule): Molecule {
    return molecule.removeAtom(this.atom.id)
  }
  description = `Add Atom ${this.atom.element}`
}

export class RemoveAtomCommand implements Command {
  private atomId: string
  private atom: Atom | null = null
  constructor(atomId: string) {
    this.atomId = atomId
  }
  execute(molecule: Molecule): Molecule {
    const atom = molecule.getAtom(this.atomId)
    if (atom) {
      this.atom = atom.toJSON()
    }
    return molecule.removeAtom(this.atomId)
  }
  undo(molecule: Molecule): Molecule {
    if (this.atom) {
      return molecule.addAtom(this.atom)
    }
    return molecule
  }
  description = `Remove Atom`
}

export class AddBondCommand implements Command {
  private bondId: string
  private atom1Id: string
  private atom2Id: string
  private order: number
  constructor(bondId: string, atom1Id: string, atom2Id: string, order: number) {
    this.bondId = bondId
    this.atom1Id = atom1Id
    this.atom2Id = atom2Id
    this.order = order
  }
  execute(molecule: Molecule): Molecule {
    return molecule.addBond({
      id: this.bondId,
      atom1: this.atom1Id,
      atom2: this.atom2Id,
      order: this.order,
    })
  }
  undo(molecule: Molecule): Molecule {
    return molecule.removeBond(this.bondId)
  }
  description = `Add Bond`
}

export class RemoveBondCommand implements Command {
  private bondId: string
  private bond: Bond | null = null
  constructor(bondId: string) {
    this.bondId = bondId
  }
  execute(molecule: Molecule): Molecule {
    const bond = molecule.getBond(this.bondId)
    if (bond) {
      this.bond = bond.toJSON()
    }
    return molecule.removeBond(this.bondId)
  }
  undo(molecule: Molecule): Molecule {
    if (this.bond) {
      return molecule.addBond(this.bond)
    }
    return molecule
  }
  description = `Remove Bond`
}

export class MoveAtomCommand implements Command {
  private atomId: string
  private oldPosition: [number, number, number]
  private newPosition: [number, number, number]
  constructor(atomId: string, oldPosition: [number, number, number], newPosition: [number, number, number]) {
    this.atomId = atomId
    this.oldPosition = oldPosition
    this.newPosition = newPosition
  }
  execute(molecule: Molecule): Molecule {
    return molecule.updateAtom(this.atomId, { position: this.newPosition })
  }
  undo(molecule: Molecule): Molecule {
    return molecule.updateAtom(this.atomId, { position: this.oldPosition })
  }
  description = `Move Atom`
}

export class UpdateAtomCommand implements Command {
  private atomId: string
  private updates: Partial<Atom>
  private previousAtom: Atom | null = null
  constructor(atomId: string, updates: Partial<Atom>) {
    this.atomId = atomId
    this.updates = updates
  }
  execute(molecule: Molecule): Molecule {
    const atom = molecule.getAtom(this.atomId)
    if (atom) {
      this.previousAtom = atom.toJSON()
    }
    return molecule.updateAtom(this.atomId, this.updates)
  }
  undo(molecule: Molecule): Molecule {
    if (this.previousAtom) {
      return molecule.updateAtom(this.atomId, this.previousAtom)
    }
    return molecule
  }
  description = `Update Atom`
}

export class UpdateBondCommand implements Command {
  private bondId: string
  private updates: Partial<Bond>
  private previousBond: Bond | null = null
  constructor(bondId: string, updates: Partial<Bond>) {
    this.bondId = bondId
    this.updates = updates
  }
  execute(molecule: Molecule): Molecule {
    const bond = molecule.getBond(this.bondId)
    if (bond) {
      this.previousBond = bond.toJSON()
    }
    return molecule.updateBond(this.bondId, this.updates)
  }
  undo(molecule: Molecule): Molecule {
    if (this.previousBond) {
      return molecule.updateBond(this.bondId, this.previousBond)
    }
    return molecule
  }
  description = `Update Bond`
}

export class ClearMoleculeCommand implements Command {
  private previousState: { atoms: Atom[]; bonds: Bond[] } | null = null
  execute(molecule: Molecule): Molecule {
    // Save current state
    this.previousState = {
      atoms: molecule.getAtoms().map(a => a.toJSON()),
      bonds: molecule.getBonds().map(b => b.toJSON()),
    }
    // Clear molecule
    const { Molecule } = require('../Molecule')
    return new Molecule()
  }
  undo(molecule: Molecule): Molecule {
    if (!this.previousState) return molecule
    const { Molecule } = require('../Molecule')
    const restored = new Molecule()
    this.previousState.atoms.forEach(atom => restored.addAtom(atom))
    this.previousState.bonds.forEach(bond => restored.addBond(bond))
    return restored
  }
  description = `Clear Molecule`
}
