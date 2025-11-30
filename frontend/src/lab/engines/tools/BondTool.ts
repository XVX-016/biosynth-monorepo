/**
 * BondTool - Tool for creating bonds between atoms
 */

import { ITool, ToolContext } from './Tool.types';

export function BondTool(stateEngine: any): ITool {
  let firstAtom: string | null = null;

  return {
    id: 'bond',
    label: 'Bond',
    icon: '🔗',

    onPointerDown(evt, ctx: ToolContext) {
      const hit = ctx.pickAtom(evt);
      if (!hit) {
        // Reset if clicking on empty space
        firstAtom = null;
        return;
      }

      if (!firstAtom) {
        // First click - select atom
        firstAtom = hit;
        // TODO: Visual feedback for selected atom
      } else {
        // Second click - create bond
        if (firstAtom !== hit) {
          try {
            stateEngine.dispatch('addBond', {
              atom1: firstAtom,
              atom2: hit,
              order: 1,
            });
          } catch (error) {
            console.warn('Failed to create bond:', error);
          }
        }
        firstAtom = null;
      }
    },

    onDeselect() {
      firstAtom = null;
    },
  };
}

