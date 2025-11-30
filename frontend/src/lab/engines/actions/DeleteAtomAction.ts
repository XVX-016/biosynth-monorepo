/**
 * DeleteAtomAction - Action for deleting an atom and its bonds
 */

import { IAction } from './Action.types';

export function DeleteAtomAction(payload: {
  atomId: string;
  atom: any;
  connectedBonds: any[];
}): IAction {
  return {
    id: crypto.randomUUID(),
    timestamp: Date.now(),
    label: `Delete Atom ${payload.atom.element}`,
    do: (state) => {
      const newAtoms = new Map(state.atoms);
      const newBonds = new Map(state.bonds);
      
      // Remove atom
      newAtoms.delete(payload.atomId);
      
      // Remove connected bonds
      payload.connectedBonds.forEach(bond => {
        newBonds.delete(bond.id);
      });
      
      return {
        ...state,
        atoms: newAtoms,
        bonds: newBonds,
      };
    },
    undo: (state) => {
      const newAtoms = new Map(state.atoms);
      const newBonds = new Map(state.bonds);
      
      // Restore atom
      newAtoms.set(payload.atomId, payload.atom);
      
      // Restore bonds
      payload.connectedBonds.forEach(bond => {
        newBonds.set(bond.id, bond);
      });
      
      return {
        ...state,
        atoms: newAtoms,
        bonds: newBonds,
      };
    },
    toJSON: () => ({
      type: 'deleteAtom',
      ...payload,
    }),
  };
}

