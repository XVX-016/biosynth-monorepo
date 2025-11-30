/**
 * EraserTool - Tool for deleting atoms and bonds
 */

import { ITool, ToolContext } from './Tool.types';

export function EraserTool(stateEngine: any): ITool {
  return {
    id: 'eraser',
    label: 'Erase',
    icon: '🗑️',

    onPointerDown(evt, ctx: ToolContext) {
      // Try to pick atom first (atoms have priority)
      const atom = ctx.pickAtom(evt);
      const bond = ctx.pickBond(evt);

      if (atom) {
        // Delete atom (will also delete connected bonds)
        const atomData = stateEngine.getAtom(atom);
        if (atomData) {
          const connectedBonds = stateEngine.getAtomBonds(atom);
          stateEngine.dispatch('deleteAtom', {
            atomId: atom,
            atom: atomData,
            connectedBonds: connectedBonds,
          });
        }
      } else if (bond) {
        // Delete bond only
        const bondData = stateEngine.getBond(bond);
        if (bondData) {
          stateEngine.dispatch('deleteBond', {
            bondId: bond,
            bond: bondData,
          });
        }
      }
    },
  };
}

