/**
 * Element data constants
 * 
 * Standard element information for validation and display.
 */

import type { ElementInfo } from './types'

export const ELEMENT_DATA: Record<string, ElementInfo> = {
  H: { symbol: 'H', atomicNumber: 1, name: 'Hydrogen', maxValence: 1, commonValence: [1] },
  He: { symbol: 'He', atomicNumber: 2, name: 'Helium', maxValence: 0, commonValence: [0] },
  Li: { symbol: 'Li', atomicNumber: 3, name: 'Lithium', maxValence: 1, commonValence: [1] },
  Be: { symbol: 'Be', atomicNumber: 4, name: 'Beryllium', maxValence: 2, commonValence: [2] },
  B: { symbol: 'B', atomicNumber: 5, name: 'Boron', maxValence: 3, commonValence: [3] },
  C: { symbol: 'C', atomicNumber: 6, name: 'Carbon', maxValence: 4, commonValence: [4] },
  N: { symbol: 'N', atomicNumber: 7, name: 'Nitrogen', maxValence: 3, commonValence: [3] },
  O: { symbol: 'O', atomicNumber: 8, name: 'Oxygen', maxValence: 2, commonValence: [2] },
  F: { symbol: 'F', atomicNumber: 9, name: 'Fluorine', maxValence: 1, commonValence: [1] },
  Ne: { symbol: 'Ne', atomicNumber: 10, name: 'Neon', maxValence: 0, commonValence: [0] },
  Na: { symbol: 'Na', atomicNumber: 11, name: 'Sodium', maxValence: 1, commonValence: [1] },
  Mg: { symbol: 'Mg', atomicNumber: 12, name: 'Magnesium', maxValence: 2, commonValence: [2] },
  Al: { symbol: 'Al', atomicNumber: 13, name: 'Aluminum', maxValence: 3, commonValence: [3] },
  Si: { symbol: 'Si', atomicNumber: 14, name: 'Silicon', maxValence: 4, commonValence: [4] },
  P: { symbol: 'P', atomicNumber: 15, name: 'Phosphorus', maxValence: 5, commonValence: [3, 5] },
  S: { symbol: 'S', atomicNumber: 16, name: 'Sulfur', maxValence: 6, commonValence: [2, 4, 6] },
  Cl: { symbol: 'Cl', atomicNumber: 17, name: 'Chlorine', maxValence: 1, commonValence: [1] },
  Ar: { symbol: 'Ar', atomicNumber: 18, name: 'Argon', maxValence: 0, commonValence: [0] },
  K: { symbol: 'K', atomicNumber: 19, name: 'Potassium', maxValence: 1, commonValence: [1] },
  Ca: { symbol: 'Ca', atomicNumber: 20, name: 'Calcium', maxValence: 2, commonValence: [2] },
  Br: { symbol: 'Br', atomicNumber: 35, name: 'Bromine', maxValence: 1, commonValence: [1] },
  I: { symbol: 'I', atomicNumber: 53, name: 'Iodine', maxValence: 1, commonValence: [1] },
}

/**
 * Common elements for the element palette
 */
export const COMMON_ELEMENTS = ['C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'H'] as const

/**
 * Default bond order
 */
export const DEFAULT_BOND_ORDER = 1

/**
 * Valid bond orders
 */
export const VALID_BOND_ORDERS = [1, 2, 3, 1.5] as const

