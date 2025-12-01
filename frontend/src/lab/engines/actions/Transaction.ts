/**
 * Transaction - Atomic multi-step action grouping
 * 
 * Ensures all actions in a transaction are applied atomically
 * and can be undone/redone as a single unit.
 */

import { IAction } from './Action.types';

export interface TransactionOptions {
  label?: string;
  validateBeforeCommit?: boolean;
  rollbackOnError?: boolean;
}

export class Transaction implements IAction {
  id: string;
  timestamp: number;
  label: string;
  private actions: IAction[] = [];
  private committed: boolean = false;
  private options: Required<TransactionOptions>;

  constructor(
    actions: IAction[] = [],
    options: TransactionOptions = {}
  ) {
    this.id = crypto.randomUUID();
    this.timestamp = Date.now();
    this.actions = [...actions];
    this.options = {
      label: options.label || `Transaction (${actions.length} actions)`,
      validateBeforeCommit: options.validateBeforeCommit ?? true,
      rollbackOnError: options.rollbackOnError ?? true,
    };
    this.label = this.options.label;
  }

  /**
   * Add an action to the transaction
   */
  addAction(action: IAction): void {
    if (this.committed) {
      throw new Error('Cannot add actions to committed transaction');
    }
    this.actions.push(action);
  }

  /**
   * Execute all actions in the transaction
   */
  do(state: any): any {
    if (this.committed) {
      throw new Error('Transaction already committed');
    }

    let currentState = state;
    const executedActions: IAction[] = [];

    try {
      // Execute all actions
      for (const action of this.actions) {
        currentState = action.do(currentState);
        executedActions.push(action);
      }

      this.committed = true;
      return currentState;
    } catch (error) {
      // Rollback executed actions if option enabled
      if (this.options.rollbackOnError) {
        // Undo in reverse order
        for (let i = executedActions.length - 1; i >= 0; i--) {
          try {
            currentState = executedActions[i].undo(currentState);
          } catch (rollbackError) {
            console.error('Rollback failed for action:', executedActions[i], rollbackError);
          }
        }
      }
      throw error;
    }
  }

  /**
   * Undo all actions in reverse order
   */
  undo(state: any): any {
    let currentState = state;

    // Undo in reverse order
    for (let i = this.actions.length - 1; i >= 0; i--) {
      currentState = this.actions[i].undo(currentState);
    }

    return currentState;
  }

  /**
   * Get transaction metadata
   */
  toJSON(): any {
    return {
      id: this.id,
      timestamp: this.timestamp,
      label: this.label,
      type: 'transaction',
      actionCount: this.actions.length,
      actions: this.actions.map(a => a.toJSON()),
    };
  }

  /**
   * Get all actions in this transaction
   */
  getActions(): IAction[] {
    return [...this.actions];
  }

  /**
   * Check if transaction is empty
   */
  isEmpty(): boolean {
    return this.actions.length === 0;
  }
}

