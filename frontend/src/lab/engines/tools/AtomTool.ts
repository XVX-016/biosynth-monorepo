/**
 * AtomTool - Tool for placing atoms
 */

import { ITool, ToolContext } from './Tool.types';

export function AtomTool(stateEngine: any): ITool {
  return {
    id: 'atom',
    label: 'Atom',
    icon: '⚛️',

    onPointerDown(evt, ctx: ToolContext) {
      const { x, y } = ctx.toMoleculeCoords(evt);
      const element = ctx.currentElement || 'C';
      
      const atomId = crypto.randomUUID();
      
      stateEngine.dispatch('addAtom', {
        atomId,
        element,
        position: { x, y },
      });
    },
  };
}

