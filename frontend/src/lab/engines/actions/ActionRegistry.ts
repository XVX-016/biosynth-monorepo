/**
 * ActionRegistry - Registry for action type builders
 */

import { IAction, ActionPayload } from './Action.types';

export class ActionRegistry {
  private registry: Record<string, (payload: ActionPayload) => IAction> = {};

  /**
   * Register an action type builder
   */
  register(type: string, builder: (payload: ActionPayload) => IAction): void {
    this.registry[type] = builder;
  }

  /**
   * Create an action instance from type and payload
   */
  create(type: string, payload: ActionPayload): IAction {
    const builder = this.registry[type];
    if (!builder) {
      throw new Error(`Unregistered action type: ${type}`);
    }
    return builder(payload);
  }

  /**
   * Check if an action type is registered
   */
  has(type: string): boolean {
    return type in this.registry;
  }

  /**
   * Get all registered action types
   */
  getTypes(): string[] {
    return Object.keys(this.registry);
  }
}

export const actionRegistry = new ActionRegistry();

