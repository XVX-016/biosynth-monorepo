/**
 * valence.ts
 * Simple valence table and helpers for suggested bond order and proximity check.
 *
 * Note: This is intentionally conservative and deterministic.
 * Chemist review recommended to extend element list and pair-wise thresholds.
 */

import type { Atom } from '@biosynth/engine';

export const STANDARD_VALENCE: Record<string, number> = {
  H: 1,
  Li: 1,
  Be: 2,
  B: 3,
  C: 4,
  N: 3,
  O: 2,
  F: 1,
  P: 3,
  S: 2,
  Cl: 1,
  Br: 1,
  I: 1,
};

export function getValence(element: string): number {
  return STANDARD_VALENCE[element] ?? 0;
}

export function distance(a: Atom, b: Atom): number {
  const dx = a.position[0] - b.position[0];
  const dy = a.position[1] - b.position[1];
  const dz = a.position[2] - b.position[2];
  return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

/**
 * Default single-bond thresholds (Å) per pair (symmetric).
 * If not defined, fallback to 1.6Å
 */
export const DEFAULT_BOND_DISTANCE = 1.6;

/**
 * Calculate current bond count for an atom
 */
export function getCurrentBondCount(atom: Atom, bonds: Array<{ a1: string; a2: string; order: number }>): number {
  return bonds.reduce((count, bond) => {
    if (bond.a1 === atom.id || bond.a2 === atom.id) {
      return count + bond.order;
    }
    return count;
  }, 0);
}

/**
 * Suggest bond order based on remaining valence
 * Conservative approach: defaults to single bond unless both atoms clearly need higher order
 */
export function suggestBondOrder(
  a: Atom,
  b: Atom,
  bonds: Array<{ a1: string; a2: string; order: number }>
): 1 | 2 | 3 {
  const currentBondsA = getCurrentBondCount(a, bonds);
  const currentBondsB = getCurrentBondCount(b, bonds);
  const remA = Math.max(0, getValence(a.element) - currentBondsA);
  const remB = Math.max(0, getValence(b.element) - currentBondsB);
  const minRem = Math.min(remA, remB);
  
  // Conservative: prefer single bonds unless there's strong evidence for higher order
  // Only suggest double/triple if both atoms have limited remaining valence
  if (minRem >= 3 && remA <= 3 && remB <= 3) return 3;
  if (minRem >= 2 && remA <= 2 && remB <= 2) return 2;
  // Default to single bond
  return 1;
}

/**
 * decideIfShouldFormBond determines whether atoms can bond given distance threshold and valence.
 */
export function shouldFormBond(
  a: Atom,
  b: Atom,
  bonds: Array<{ a1: string; a2: string; order: number }>,
  threshold = DEFAULT_BOND_DISTANCE
): boolean {
  const d = distance(a, b);
  if (d > threshold) return false;
  const currentBondsA = getCurrentBondCount(a, bonds);
  const currentBondsB = getCurrentBondCount(b, bonds);
  const remA = Math.max(0, getValence(a.element) - currentBondsA);
  const remB = Math.max(0, getValence(b.element) - currentBondsB);
  return remA >= 1 && remB >= 1;
}

