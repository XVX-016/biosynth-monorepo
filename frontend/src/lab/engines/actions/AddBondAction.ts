/**
 * AddBondAction - Action for adding a bond
 */

import { IAction, ActionPayload } from './Action.types';

export function AddBondAction(payload: {
  bondId: string;
  atom1: string;
  atom2: string;
  order: number;
}): IAction {
  return {
    id: crypto.randomUUID(),
    timestamp: Date.now(),
    label: `Add Bond (order ${payload.order})`,
    do: (state) => {
      const newBonds = new Map(state.bonds);
      newBonds.set(payload.bondId, {
        id: payload.bondId,
        atoms: [payload.atom1, payload.atom2],
        order: payload.order,
      });
      return {
        ...state,
        bonds: newBonds,
      };
    },
    undo: (state) => {
      const newBonds = new Map(state.bonds);
      newBonds.delete(payload.bondId);
      return {
        ...state,
        bonds: newBonds,
      };
    },
    toJSON: () => ({
      type: 'addBond',
      ...payload,
    }),
  };
}

