/**
 * AddAtomAction - Action for adding an atom
 */

import { IAction, ActionPayload } from './Action.types';

export function AddAtomAction(payload: {
  atomId: string;
  element: string;
  position: { x: number; y: number; z?: number };
  charge?: number;
}): IAction {
  return {
    id: crypto.randomUUID(),
    timestamp: Date.now(),
    label: `Add Atom ${payload.element}`,
    do: (state) => {
      const newAtoms = new Map(state.atoms);
      newAtoms.set(payload.atomId, {
        id: payload.atomId,
        element: payload.element,
        x: payload.position.x,
        y: payload.position.y,
        z: payload.position.z || 0,
        charge: payload.charge || 0,
      });
      return {
        ...state,
        atoms: newAtoms,
      };
    },
    undo: (state) => {
      const newAtoms = new Map(state.atoms);
      newAtoms.delete(payload.atomId);
      return {
        ...state,
        atoms: newAtoms,
      };
    },
    toJSON: () => ({
      type: 'addAtom',
      ...payload,
    }),
  };
}

