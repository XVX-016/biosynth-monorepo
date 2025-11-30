/**
 * Action System - Index
 * 
 * Registers all action types with the registry
 */

import { actionRegistry } from './ActionRegistry';
import { AddAtomAction } from './AddAtomAction';
import { AddBondAction } from './AddBondAction';
import { DeleteAtomAction } from './DeleteAtomAction';
import { DeleteBondAction } from './DeleteBondAction';
import { MoveAtomAction } from './MoveAtomAction';

/**
 * Register all action types
 */
export function registerActions(): void {
  actionRegistry.register('addAtom', AddAtomAction);
  actionRegistry.register('addBond', AddBondAction);
  actionRegistry.register('deleteAtom', DeleteAtomAction);
  actionRegistry.register('deleteBond', DeleteBondAction);
  actionRegistry.register('moveAtom', MoveAtomAction);
}

// Auto-register on import
registerActions();

// Re-export everything
export * from './Action.types';
export * from './ActionManager';
export * from './ActionRegistry';
export * from './AddAtomAction';
export * from './AddBondAction';
export * from './DeleteAtomAction';
export * from './DeleteBondAction';
export * from './MoveAtomAction';

