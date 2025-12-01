/**
 * ActionManager - Manages undo/redo stacks and action execution
 * 
 * Enhanced with transaction support for atomic multi-step operations.
 */

import { IAction } from './Action.types';
import { Transaction } from './Transaction';

export class ActionManager {
  private undoStack: IAction[] = [];
  private redoStack: IAction[] = [];
  private listeners: ((state: any) => void)[] = [];
  private maxStackSize = 100;
  private currentTransaction: Transaction | null = null;
  private transactionDepth = 0;

  /**
   * Begin a transaction for atomic multi-step operations
   */
  beginTransaction(label?: string): void {
    if (this.currentTransaction === null) {
      this.currentTransaction = new Transaction([], { label });
      this.transactionDepth = 1;
    } else {
      // Nested transaction - will be added to parent
      this.transactionDepth++;
    }
  }

  /**
   * Commit the current transaction
   */
  commitTransaction(state: any): any {
    if (this.currentTransaction === null) {
      throw new Error('No active transaction to commit');
    }

    this.transactionDepth--;
    
    if (this.transactionDepth > 0) {
      // Nested transaction - return state unchanged, parent will handle
      return state;
    }

    // Commit the transaction
    const transaction = this.currentTransaction;
    this.currentTransaction = null;

    if (transaction.isEmpty()) {
      return state; // Empty transaction, no-op
    }

    // Apply transaction as a single action
    return this.apply(transaction, state);
  }

  /**
   * Rollback the current transaction
   */
  rollbackTransaction(state: any): any {
    if (this.currentTransaction === null) {
      throw new Error('No active transaction to rollback');
    }

    this.transactionDepth--;
    
    if (this.transactionDepth > 0) {
      // Nested transaction - clear current, let parent handle
      this.currentTransaction = null;
      return state;
    }

    // Rollback all actions in transaction
    const transaction = this.currentTransaction;
    this.currentTransaction = null;

    if (transaction.isEmpty()) {
      return state; // Empty transaction, no-op
    }

    // Undo all actions in reverse order
    let currentState = state;
    const actions = transaction.getActions();
    for (let i = actions.length - 1; i >= 0; i--) {
      currentState = actions[i].undo(currentState);
    }

    return currentState;
  }

  /**
   * Check if a transaction is active
   */
  isInTransaction(): boolean {
    return this.currentTransaction !== null;
  }

  /**
   * Apply an action and add it to undo stack
   * If in a transaction, adds to transaction instead
   */
  apply(action: IAction, state: any): any {
    // If in transaction, add to transaction instead of applying immediately
    if (this.currentTransaction !== null) {
      this.currentTransaction.addAction(action);
      // Execute action immediately but don't add to undo stack yet
      const newState = action.do(state);
      return newState;
    }

    // Normal application
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

