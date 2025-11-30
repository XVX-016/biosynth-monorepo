/**
 * BasicValenceValidator - Validates atom valences
 */

import { ValidationError } from './Validation.types';
import type { MoleculeState } from '../MoleculeStateEngine';

const VALENCE_TABLE: Record<string, number> = {
  H: 1,
  C: 4,
  N: 3,
  O: 2,
  F: 1,
  S: 6,
  Cl: 1,
  Br: 1,
  I: 1,
  P: 5,
  B: 3,
  Si: 4,
};

export function validateValence(state: MoleculeState): ValidationError[] {
  const errors: ValidationError[] = [];

  state.atoms.forEach((atom, atomId) => {
    // Get all bonds connected to this atom
    const connectedBonds = Array.from(state.bonds.values()).filter(
      (bond) => bond.atoms[0] === atomId || bond.atoms[1] === atomId
    );

    // Sum up bond orders
    const bondOrderSum = connectedBonds.reduce(
      (acc, bond) => acc + (bond.order || 1),
      0
    );

    // Get expected valence
    const expectedValence = VALENCE_TABLE[atom.element] ?? 4;

    // Check if exceeds valence
    if (bondOrderSum > expectedValence) {
      errors.push({
        type: 'valence',
        message: `Atom ${atom.element} exceeds allowed valence (${bondOrderSum} > ${expectedValence})`,
        atomId,
      });
    }

    // Check if under-bonded (warning, not error)
    if (bondOrderSum < expectedValence && atom.element !== 'H') {
      // This is just a warning, not an error
      // Some atoms can be under-bonded (e.g., radicals)
    }
  });

  return errors;
}

