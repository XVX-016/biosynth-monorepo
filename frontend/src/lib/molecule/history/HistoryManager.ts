/**
 * HistoryManager - Manages undo/redo history
 * 
 * Phase 4: Undo/Redo System
 * 
 * Maintains a command stack with undo/redo capabilities.
 */

import type { Molecule } from '../Molecule'
import type { Command } from './Command'

export class HistoryManager {
  private history: Command[] = []
  private future: Command[] = []
  private maxHistorySize: number = 100

  /**
   * Execute a command and add it to history
   */
  execute(molecule: Molecule, command: Command): void {
    command.execute(molecule)
    this.history.push(command)
    
    // Clear future when new command is executed
    this.future = []
    
    // Limit history size
    if (this.history.length > this.maxHistorySize) {
      this.history.shift()
    }
  }

  /**
   * Undo last command
   */
  undo(molecule: Molecule): boolean {
    if (this.history.length === 0) {
      return false
    }

    const command = this.history.pop()!
    command.undo(molecule)
    this.future.push(command)
    return true
  }

  /**
   * Redo last undone command
   */
  redo(molecule: Molecule): boolean {
    if (this.future.length === 0) {
      return false
    }

    const command = this.future.pop()!
    command.execute(molecule)
    this.history.push(command)
    return true
  }

  /**
   * Check if undo is available
   */
  canUndo(): boolean {
    return this.history.length > 0
  }

  /**
   * Check if redo is available
   */
  canRedo(): boolean {
    return this.future.length > 0
  }

  /**
   * Get undo description
   */
  getUndoDescription(): string | null {
    if (this.history.length === 0) return null
    return this.history[this.history.length - 1].getDescription()
  }

  /**
   * Get redo description
   */
  getRedoDescription(): string | null {
    if (this.future.length === 0) return null
    return this.future[this.future.length - 1].getDescription()
  }

  /**
   * Clear all history
   */
  clear(): void {
    this.history = []
    this.future = []
  }

  /**
   * Get history size
   */
  getHistorySize(): number {
    return this.history.length
  }

  /**
   * Get future size
   */
  getFutureSize(): number {
    return this.future.length
  }

  /**
   * Set maximum history size
   */
  setMaxHistorySize(size: number): void {
    this.maxHistorySize = size
    // Trim history if needed
    while (this.history.length > this.maxHistorySize) {
      this.history.shift()
    }
  }
}

