/**
 * SelectTool - Tool for selecting atoms and bonds
 */

import { ITool, ToolContext } from './Tool.types';

export function SelectTool(stateEngine: any, selectionManager: any): ITool {
  return {
    id: 'select',
    label: 'Select',
    icon: '👆',

    onPointerDown(evt, ctx: ToolContext) {
      const atom = ctx.pickAtom(evt);
      const bond = ctx.pickBond(evt);

      if (atom) {
        selectionManager.selectAtom(atom);
      } else if (bond) {
        selectionManager.selectBond(bond);
      } else {
        selectionManager.clear();
      }
    },
  };
}

