/**
 * bonds.ts
 * Kernel-level operations to manage bonds with validation.
 * Works with @biosynth/engine MoleculeGraph structure.
 */

import { MoleculeGraph } from '@biosynth/engine';
import { suggestBondOrder, shouldFormBond } from './valence';

/**
 * createBondSafe - create a bond only if valence and other checks pass.
 * Returns bond ID or null if invalid.
 */
export function createBondSafe(
  molecule: MoleculeGraph,
  aId: string,
  bId: string,
  order?: 1 | 2 | 3
): string | null {
  const a = molecule.atoms.get(aId);
  const b = molecule.atoms.get(bId);
  if (!a || !b) return null;

  // If already a bond exists, return existing bond ID
  const existingBond = Array.from(molecule.bonds.values()).find(
    (bond) => (bond.a1 === aId && bond.a2 === bId) || (bond.a1 === bId && bond.a2 === aId)
  );
  if (existingBond) return existingBond.id;

  // Get all bonds as array for valence checking
  const bondsArray = Array.from(molecule.bonds.values());

  // Use suggested order if not provided
  const suggested = order ?? suggestBondOrder(a, b, bondsArray);

  // Check simple valence/distance test
  if (!shouldFormBond(a, b, bondsArray)) {
    // If they're out of proximity but user forces via UI, KERNEL can still create (but here we block)
    return null;
  }

  // Create bond in molecule
  const bondId = molecule.addBond(aId, bId, suggested);
  return bondId;
}

/**
 * removeBondSafe - remove bond from molecule
 */
export function removeBondSafe(molecule: MoleculeGraph, bondId: string): boolean {
  return molecule.removeBond(bondId);
}

/**
 * Check if two atoms have a bond between them
 */
export function hasBondBetween(molecule: MoleculeGraph, aId: string, bId: string): boolean {
  return Array.from(molecule.bonds.values()).some(
    (bond) => (bond.a1 === aId && bond.a2 === bId) || (bond.a1 === bId && bond.a2 === aId)
  );
}

/**
 * Find bond between two atoms
 */
export function findBondBetween(
  molecule: MoleculeGraph,
  aId: string,
  bId: string
): { id: string; a1: string; a2: string; order: number } | undefined {
  return Array.from(molecule.bonds.values()).find(
    (bond) => (bond.a1 === aId && bond.a2 === bId) || (bond.a1 === bId && bond.a2 === aId)
  );
}

