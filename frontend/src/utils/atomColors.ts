/**
 * Shared Atom Color System
 * 
 * Single source of truth for atom colors across Library and Lab.
 * Uses CPK (Corey-Pauling-Koltun) color scheme for consistency.
 */

export type ElementSymbol = 'H' | 'C' | 'N' | 'O' | 'S' | 'F' | 'Cl' | 'Br' | 'P' | 'I';

/**
 * CPK colors as hex strings (for CSS/Three.js string colors)
 */
export const ATOM_COLORS: Record<string, string> = {
  H: '#FFFFFF',   // White
  C: '#909090',   // Dark gray
  N: '#3050F8',   // Blue
  O: '#FF0D0D',   // Red
  F: '#90E050',   // Green
  P: '#FF8000',   // Orange
  S: '#FFFF30',   // Yellow
  Cl: '#1FF01F',  // Green
  Br: '#A62929',  // Dark red
  I: '#940094',   // Purple
};

/**
 * CPK colors as numbers (for Three.js numeric colors)
 */
export const ATOM_COLORS_NUMERIC: Record<string, number> = {
  H: 0xffffff,
  C: 0x909090,
  N: 0x3050f8,
  O: 0xff0d0d,
  F: 0x90e050,
  P: 0xff8000,
  S: 0xffff30,
  Cl: 0x1ff01f,
  Br: 0xa62929,
  I: 0x940094,
};

/**
 * Default color for unknown elements
 */
export const DEFAULT_ATOM_COLOR = '#888888';
export const DEFAULT_ATOM_COLOR_NUMERIC = 0x888888;

/**
 * Get atom color as hex string
 */
export function getAtomColor(element: string): string {
  return ATOM_COLORS[element] || DEFAULT_ATOM_COLOR;
}

/**
 * Get atom color as numeric (for Three.js)
 */
export function getAtomColorNumeric(element: string): number {
  return ATOM_COLORS_NUMERIC[element] || DEFAULT_ATOM_COLOR_NUMERIC;
}

/**
 * Get text color for atom (for contrast on colored backgrounds)
 */
export function getAtomTextColor(element: string): string {
  // Yellow (S) needs dark text, white (H) needs dark text
  if (element === 'S' || element === 'H') {
    return '#000000';
  }
  return '#FFFFFF';
}

