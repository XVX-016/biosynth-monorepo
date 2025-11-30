/**
 * CommandStack - Implements undo/redo pattern for molecule operations
 */

import { moleculeEngine, type Atom, type Bond } from './MoleculeStateEngine'

export interface Command {
  execute(): void
  undo(): void
  description?: string
}

/**
 * Command to add an atom
 */
export class AddAtomCommand implements Command {
  private atomId: string
  private element: string
  private x: number
  private y: number
  private z: number

  constructor(element: string, x: number, y: number, z = 0) {
    this.element = element
    this.x = x
    this.y = y
    this.z = z
    this.atomId = ""
  }

  execute(): void {
    this.atomId = moleculeEngine.addAtom(this.element, this.x, this.y, this.z)
  }

  undo(): void {
    if (this.atomId) {
      moleculeEngine.removeAtom(this.atomId)
    }
  }

  get description(): string {
    return `Add ${this.element} atom`
  }
}

/**
 * Command to remove an atom
 */
export class RemoveAtomCommand implements Command {
  private atomId: string
  private atom: Atom | undefined
  private connectedBonds: Bond[] = []

  constructor(atomId: string) {
    this.atomId = atomId
  }

  execute(): void {
    this.atom = moleculeEngine.getAtom(this.atomId)
    if (this.atom) {
      // Save connected bonds for undo
      this.connectedBonds = moleculeEngine.getBondsForAtom(this.atomId)
      moleculeEngine.removeAtom(this.atomId)
    }
  }

  undo(): void {
    if (this.atom) {
      // Restore atom
      moleculeEngine.atoms.set(this.atom.id, this.atom)
      // Restore bonds
      this.connectedBonds.forEach(bond => {
        moleculeEngine.bonds.set(bond.id, bond)
      })
    }
  }

  get description(): string {
    return `Remove atom`
  }
}

/**
 * Command to add a bond
 */
export class AddBondCommand implements Command {
  private bondId: string | null = null
  private atomA: string
  private atomB: string
  private order: number

  constructor(atomA: string, atomB: string, order = 1) {
    this.atomA = atomA
    this.atomB = atomB
    this.order = order
  }

  execute(): void {
    this.bondId = moleculeEngine.addBond(this.atomA, this.atomB, this.order)
  }

  undo(): void {
    if (this.bondId) {
      moleculeEngine.removeBond(this.bondId)
    }
  }

  get description(): string {
    return `Add bond`
  }
}

/**
 * Command to move an atom
 */
export class MoveAtomCommand implements Command {
  private atomId: string
  private oldX: number
  private oldY: number
  private oldZ: number | undefined
  private newX: number
  private newY: number
  private newZ: number | undefined

  constructor(atomId: string, newX: number, newY: number, newZ?: number) {
    this.atomId = atomId
    const atom = moleculeEngine.getAtom(atomId)
    if (atom) {
      this.oldX = atom.x
      this.oldY = atom.y
      this.oldZ = atom.z
    } else {
      this.oldX = 0
      this.oldY = 0
      this.oldZ = 0
    }
    this.newX = newX
    this.newY = newY
    this.newZ = newZ
  }

  execute(): void {
    moleculeEngine.setAtomPosition(this.atomId, this.newX, this.newY, this.newZ)
  }

  undo(): void {
    moleculeEngine.setAtomPosition(this.atomId, this.oldX, this.oldY, this.oldZ)
  }

  get description(): string {
    return `Move atom`
  }
}

/**
 * Command stack manager
 */
export class CommandStack {
  private undoStack: Command[] = []
  private redoStack: Command[] = []
  private maxStackSize = 100

  /**
   * Execute a command and add to undo stack
   */
  run(cmd: Command): void {
    cmd.execute()
    this.undoStack.push(cmd)
    
    // Limit stack size
    if (this.undoStack.length > this.maxStackSize) {
      this.undoStack.shift()
    }
    
    // Clear redo stack when new command is executed
    this.redoStack = []
  }

  /**
   * Undo the last command
   */
  undo(): boolean {
    const cmd = this.undoStack.pop()
    if (cmd) {
      cmd.undo()
      this.redoStack.push(cmd)
      return true
    }
    return false
  }

  /**
   * Redo the last undone command
   */
  redo(): boolean {
    const cmd = this.redoStack.pop()
    if (cmd) {
      cmd.execute()
      this.undoStack.push(cmd)
      return true
    }
    return false
  }

  /**
   * Check if undo is available
   */
  canUndo(): boolean {
    return this.undoStack.length > 0
  }

  /**
   * Check if redo is available
   */
  canRedo(): boolean {
    return this.redoStack.length > 0
  }

  /**
   * Clear all stacks
   */
  clear(): void {
    this.undoStack = []
    this.redoStack = []
  }

  /**
   * Get undo stack size
   */
  getUndoSize(): number {
    return this.undoStack.length
  }

  /**
   * Get redo stack size
   */
  getRedoSize(): number {
    return this.redoStack.length
  }
}

// Singleton instance
export const commandStack = new CommandStack()

