/**
 * DeleteBondAction - Action for deleting a bond
 */

import { IAction } from './Action.types';

export function DeleteBondAction(payload: {
  bondId: string;
  bond: any;
}): IAction {
  return {
    id: crypto.randomUUID(),
    timestamp: Date.now(),
    label: 'Delete Bond',
    do: (state) => {
      const newBonds = new Map(state.bonds);
      newBonds.delete(payload.bondId);
      return {
        ...state,
        bonds: newBonds,
      };
    },
    undo: (state) => {
      const newBonds = new Map(state.bonds);
      newBonds.set(payload.bondId, payload.bond);
      return {
        ...state,
        bonds: newBonds,
      };
    },
    toJSON: () => ({
      type: 'deleteBond',
      ...payload,
    }),
  };
}

