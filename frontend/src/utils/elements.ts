/**
 * Element Specification System
 * 
 * Defines valency rules and properties for supported elements.
 * Pragmatic approach: supports common organic chemistry elements only.
 */

export type ElementSymbol = 'H' | 'C' | 'N' | 'O' | 'S' | 'F' | 'Cl' | 'Br';

export type BondOrder = 1 | 2 | 3;

export interface ElementSpec {
  symbol: ElementSymbol;
  allowedValence: number[];
  maxBondOrder: BondOrder;
  radius: number; // Van der Waals radius in Angstroms (for rendering)
}

/**
 * Element specifications with valency rules
 * 
 * Notes:
 * - N: Can have valence 3 (amines) or 5 (ammonium, nitro groups)
 * - S: Simplified to 2/4/6 (common oxidation states)
 * - Halogens: Typically single bonds only
 */
export const ELEMENTS: Record<ElementSymbol, ElementSpec> = {
  H: {
    symbol: 'H',
    allowedValence: [1],
    maxBondOrder: 1,
    radius: 0.37,
  },
  C: {
    symbol: 'C',
    allowedValence: [4],
    maxBondOrder: 3,
    radius: 0.77,
  },
  N: {
    symbol: 'N',
    allowedValence: [3, 5],
    maxBondOrder: 3,
    radius: 0.75,
  },
  O: {
    symbol: 'O',
    allowedValence: [2],
    maxBondOrder: 2,
    radius: 0.73,
  },
  S: {
    symbol: 'S',
    allowedValence: [2, 4, 6],
    maxBondOrder: 2, // Simplified: S=S double bonds, but not S≡S triple
    radius: 1.02,
  },
  F: {
    symbol: 'F',
    allowedValence: [1],
    maxBondOrder: 1,
    radius: 0.71,
  },
  Cl: {
    symbol: 'Cl',
    allowedValence: [1],
    maxBondOrder: 1,
    radius: 0.99,
  },
  Br: {
    symbol: 'Br',
    allowedValence: [1],
    maxBondOrder: 1,
    radius: 1.14,
  },
};

/**
 * Get element specification
 */
export function getElementSpec(element: string): ElementSpec | null {
  return ELEMENTS[element as ElementSymbol] || null;
}

/**
 * Check if an element is supported
 */
export function isElementSupported(element: string): boolean {
  return element in ELEMENTS;
}

/**
 * Get all supported element symbols
 */
export function getSupportedElements(): ElementSymbol[] {
  return Object.keys(ELEMENTS) as ElementSymbol[];
}

