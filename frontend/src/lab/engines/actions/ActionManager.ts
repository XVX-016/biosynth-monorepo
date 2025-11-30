/**
 * ActionManager - Manages undo/redo stacks and action execution
 */

import { IAction } from './Action.types';

export class ActionManager {
  private undoStack: IAction[] = [];
  private redoStack: IAction[] = [];
  private listeners: ((state: any) => void)[] = [];
  private maxStackSize = 100;

  /**
   * Apply an action and add it to undo stack
   */
  apply(action: IAction, state: any): any {
    const newState = action.do(state);
    this.undoStack.push(action);
    
    // Limit stack size
    if (this.undoStack.length > this.maxStackSize) {
      this.undoStack.shift();
    }
    
    this.redoStack = [];
    this.notify(newState);
    return newState;
  }

  /**
   * Undo the last action
   */
  undo(state: any): any {
    const action = this.undoStack.pop();
    if (!action) return state;
    
    const newState = action.undo(state);
    this.redoStack.push(action);
    this.notify(newState);
    return newState;
  }

  /**
   * Redo the last undone action
   */
  redo(state: any): any {
    const action = this.redoStack.pop();
    if (!action) return state;
    
    const newState = action.do(state);
    this.undoStack.push(action);
    this.notify(newState);
    return newState;
  }

  /**
   * Check if undo is available
   */
  canUndo(): boolean {
    return this.undoStack.length > 0;
  }

  /**
   * Check if redo is available
   */
  canRedo(): boolean {
    return this.redoStack.length > 0;
  }

  /**
   * Get undo stack for history display
   */
  getUndoStack(): IAction[] {
    return [...this.undoStack];
  }

  /**
   * Get redo stack
   */
  getRedoStack(): IAction[] {
    return [...this.redoStack];
  }

  /**
   * Clear all stacks
   */
  clear(): void {
    this.undoStack = [];
    this.redoStack = [];
  }

  /**
   * Subscribe to state changes
   */
  onChange(cb: (state: any) => void): () => void {
    this.listeners.push(cb);
    return () => {
      const index = this.listeners.indexOf(cb);
      if (index > -1) {
        this.listeners.splice(index, 1);
      }
    };
  }

  /**
   * Notify all listeners of state change
   */
  private notify(state: any): void {
    this.listeners.forEach((cb) => cb(state));
  }
}

