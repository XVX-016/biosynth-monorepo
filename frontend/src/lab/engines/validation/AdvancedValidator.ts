/**
 * AdvancedValidator - Extended validation rules
 */

import { ValidationError } from './Validation.types';
import type { MoleculeState } from '../MoleculeStateEngine';

/**
 * Validate formal charges
 */
export function validateCharges(state: MoleculeState): ValidationError[] {
  const errors: ValidationError[] = [];

  state.atoms.forEach((atom, atomId) => {
    // Check for unrealistic charges
    if (Math.abs(atom.charge) > 3) {
      errors.push({
        type: 'charge',
        message: `Atom ${atom.element} has unrealistic charge: ${atom.charge}`,
        atomId,
      });
    }

    // Check for charge on hydrogen
    if (atom.element === 'H' && atom.charge !== 0) {
      errors.push({
        type: 'charge',
        message: 'Hydrogen should not have formal charge',
        atomId,
      });
    }
  });

  // Check total charge
  const totalCharge = Array.from(state.atoms.values()).reduce(
    (sum, atom) => sum + atom.charge,
    0
  );
  if (Math.abs(totalCharge) > 0 && state.atoms.size > 1) {
    errors.push({
      type: 'charge',
      message: `Molecule has non-zero total charge: ${totalCharge}`,
    });
  }

  return errors;
}

/**
 * Detect aromaticity (basic check)
 */
export function validateAromaticity(state: MoleculeState): ValidationError[] {
  const errors: ValidationError[] = [];
  const warnings: ValidationError[] = [];

  // Check for potential aromatic rings (6-membered rings with alternating bonds)
  // This is a simplified check - full aromaticity detection requires RDKit
  const rings = detectRings(state);

  rings.forEach((ring, index) => {
    if (ring.length === 6) {
      // Check if all bonds are single or double (potential aromatic)
      const hasAlternatingBonds = checkAlternatingBonds(state, ring);
      if (hasAlternatingBonds) {
        warnings.push({
          type: 'aromaticity',
          message: `Ring ${index + 1} may be aromatic (requires RDKit confirmation)`,
        });
      }
    }
  });

  return [...errors, ...warnings];
}

/**
 * Validate hybridization (basic check)
 */
export function validateHybridization(state: MoleculeState): ValidationError[] {
  const errors: ValidationError[] = [];

  state.atoms.forEach((atom, atomId) => {
    const bonds = Array.from(state.bonds.values()).filter(
      (bond) => bond.atoms[0] === atomId || bond.atoms[1] === atomId
    );

    const bondCount = bonds.length;
    const bondOrderSum = bonds.reduce((sum, bond) => sum + bond.order, 0);

    // Basic hybridization checks
    if (atom.element === 'C') {
      if (bondCount === 4 && bondOrderSum === 4) {
        // sp3 - OK
      } else if (bondCount === 3 && bondOrderSum === 4) {
        // sp2 - OK
      } else if (bondCount === 2 && bondOrderSum === 4) {
        // sp - OK
      } else if (bondCount > 4) {
        errors.push({
          type: 'hybridization',
          message: `Carbon atom has ${bondCount} bonds (max 4)`,
          atomId,
        });
      }
    }

    if (atom.element === 'N') {
      if (bondCount > 4) {
        errors.push({
          type: 'hybridization',
          message: `Nitrogen atom has ${bondCount} bonds (max 4)`,
          atomId,
        });
      }
    }
  });

  return errors;
}

/**
 * Detect rings in molecule (simplified)
 */
function detectRings(state: MoleculeState): string[][] {
  const rings: string[][] = [];
  const visited = new Set<string>();

  const dfs = (start: string, current: string, path: string[], visitedInPath: Set<string>) => {
    if (current === start && path.length > 2) {
      // Found a ring
      rings.push([...path]);
      return;
    }

    if (visitedInPath.has(current)) return;

    visitedInPath.add(current);
    path.push(current);

    const bonds = Array.from(state.bonds.values()).filter(
      (bond) => bond.atoms[0] === current || bond.atoms[1] === current
    );

    bonds.forEach((bond) => {
      const next = bond.atoms[0] === current ? bond.atoms[1] : bond.atoms[0];
      if (!visitedInPath.has(next) || (next === start && path.length > 2)) {
        dfs(start, next, [...path], new Set(visitedInPath));
      }
    });
  };

  // Try starting from each atom
  state.atoms.forEach((atom) => {
    if (!visited.has(atom.id)) {
      dfs(atom.id, atom.id, [], new Set());
      visited.add(atom.id);
    }
  });

  // Remove duplicate rings
  return rings.filter((ring, index, self) => {
    const normalized = ring.sort().join('-');
    return index === self.findIndex((r) => r.sort().join('-') === normalized);
  });
}

/**
 * Check if ring has alternating single/double bonds
 */
function checkAlternatingBonds(state: MoleculeState, ring: string[]): boolean {
  if (ring.length < 3) return false;

  // Simplified check - just verify bonds exist
  for (let i = 0; i < ring.length; i++) {
    const atom1 = ring[i];
    const atom2 = ring[(i + 1) % ring.length];

    const bond = Array.from(state.bonds.values()).find(
      (b) =>
        (b.atoms[0] === atom1 && b.atoms[1] === atom2) ||
        (b.atoms[0] === atom2 && b.atoms[1] === atom1)
    );

    if (!bond) return false;
  }

  return true;
}

/**
 * Run all advanced validations
 */
export function runAdvancedValidations(state: MoleculeState): ValidationError[] {
  const errors: ValidationError[] = [];

  errors.push(...validateCharges(state));
  errors.push(...validateHybridization(state));
  // Aromaticity returns warnings, not errors
  // errors.push(...validateAromaticity(state));

  return errors;
}

