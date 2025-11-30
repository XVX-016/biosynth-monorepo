/**
 * DragTool - Tool for moving atoms
 */

import { ITool, ToolContext } from './Tool.types';

export function DragTool(stateEngine: any): ITool {
  let draggedAtom: string | null = null;
  let startPosition: { x: number; y: number } | null = null;

  return {
    id: 'drag',
    label: 'Drag',
    icon: '✋',

    onPointerDown(evt, ctx: ToolContext) {
      const hit = ctx.pickAtom(evt);
      if (hit) {
        draggedAtom = hit;
        const atom = stateEngine.getAtom(hit);
        if (atom) {
          startPosition = { x: atom.x, y: atom.y };
        }
      }
    },

    onPointerMove(evt, ctx: ToolContext) {
      if (!draggedAtom || !startPosition) return;
      
      const { x, y } = ctx.toMoleculeCoords(evt);
      
      stateEngine.dispatch('moveAtom', {
        atomId: draggedAtom,
        oldPosition: startPosition,
        newPosition: { x, y },
      });
      
      // Update start position for smooth dragging
      startPosition = { x, y };
    },

    onPointerUp() {
      draggedAtom = null;
      startPosition = null;
    },

    onDeselect() {
      draggedAtom = null;
      startPosition = null;
    },
  };
}

