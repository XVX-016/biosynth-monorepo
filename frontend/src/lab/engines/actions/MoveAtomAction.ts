/**
 * MoveAtomAction - Action for moving an atom
 */

import { IAction } from './Action.types';

export function MoveAtomAction(payload: {
  atomId: string;
  oldPosition: { x: number; y: number; z?: number };
  newPosition: { x: number; y: number; z?: number };
}): IAction {
  return {
    id: crypto.randomUUID(),
    timestamp: Date.now(),
    label: 'Move Atom',
    do: (state) => {
      const newAtoms = new Map(state.atoms);
      const atom = newAtoms.get(payload.atomId);
      if (atom) {
        newAtoms.set(payload.atomId, {
          ...atom,
          x: payload.newPosition.x,
          y: payload.newPosition.y,
          z: payload.newPosition.z ?? atom.z ?? 0,
        });
      }
      return {
        ...state,
        atoms: newAtoms,
      };
    },
    undo: (state) => {
      const newAtoms = new Map(state.atoms);
      const atom = newAtoms.get(payload.atomId);
      if (atom) {
        newAtoms.set(payload.atomId, {
          ...atom,
          x: payload.oldPosition.x,
          y: payload.oldPosition.y,
          z: payload.oldPosition.z ?? atom.z ?? 0,
        });
      }
      return {
        ...state,
        atoms: newAtoms,
      };
    },
    toJSON: () => ({
      type: 'moveAtom',
      ...payload,
    }),
  };
}

