/**
 * Tool System - Types
 */

export interface ITool {
  id: string;
  icon?: string;
  label: string;
  
  onPointerDown?(evt: any, ctx: any): void;
  onPointerMove?(evt: any, ctx: any): void;
  onPointerUp?(evt: any, ctx: any): void;
  
  onSelect?(): void;
  onDeselect?(): void;
}

export interface ToolContext {
  toMoleculeCoords(evt: any): { x: number; y: number };
  pickAtom(evt: any): string | null;
  pickBond(evt: any): string | null;
  currentElement?: string;
  stateEngine?: any;
}

