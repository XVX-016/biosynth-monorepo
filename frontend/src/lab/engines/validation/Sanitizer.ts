/**
 * Sanitizer - Cleans up molecule state
 */

import type { MoleculeState, Atom, Bond } from '../MoleculeStateEngine';

/**
 * Normalize hydrogens - add implicit hydrogens based on valence
 */
function normalizeHydrogens(state: MoleculeState): MoleculeState {
  const newState: MoleculeState = {
    atoms: new Map(state.atoms),
    bonds: new Map(state.bonds),
  };

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

  // For each atom, check if it needs hydrogens
  newState.atoms.forEach((atom, atomId) => {
    if (atom.element === 'H') return; // Skip hydrogen atoms

    const bonds = Array.from(newState.bonds.values()).filter(
      (bond) => bond.atoms[0] === atomId || bond.atoms[1] === atomId
    );

    const bondOrderSum = bonds.reduce((sum, bond) => sum + bond.order, 0);
    const expectedValence = VALENCE_TABLE[atom.element] || 4;
    const missingValence = expectedValence - bondOrderSum;

    // Add implicit hydrogens (for display purposes, we could add them as atoms)
    // For now, we just track the information
    if (missingValence > 0 && missingValence <= expectedValence) {
      // This atom needs hydrogens, but we'll leave it for now
      // In a full implementation, we'd add H atoms here
    }
  });

  return newState;
}

/**
 * Fix formal charges - normalize charges based on bonding
 */
function fixCharges(state: MoleculeState): MoleculeState {
  const newState: MoleculeState = {
    atoms: new Map(state.atoms),
    bonds: new Map(state.bonds),
  };

  // Reset charges to 0 for now (full implementation would calculate based on bonding)
  newState.atoms.forEach((atom) => {
    if (Math.abs(atom.charge) > 3) {
      atom.charge = 0; // Reset unrealistic charges
    }
  });

  return newState;
}

/**
 * Remove duplicate bonds
 */
function removeDuplicateBonds(state: MoleculeState): MoleculeState {
  const newState: MoleculeState = {
    atoms: new Map(state.atoms),
    bonds: new Map(),
  };

  const seen = new Set<string>();

  state.bonds.forEach((bond) => {
    const key = [bond.atoms[0], bond.atoms[1]].sort().join('-');
    if (!seen.has(key)) {
      seen.add(key);
      newState.bonds.set(bond.id, bond);
    }
  });

  return newState;
}

/**
 * Clean up geometry - ensure atoms aren't too close together
 */
function cleanupGeometry(state: MoleculeState): MoleculeState {
  const newState: MoleculeState = {
    atoms: new Map(state.atoms),
    bonds: new Map(state.bonds),
  };

  const MIN_DISTANCE = 1.0; // Minimum distance between atoms

  const atoms = Array.from(newState.atoms.values());
  for (let i = 0; i < atoms.length; i++) {
    for (let j = i + 1; j < atoms.length; j++) {
      const dist = Math.sqrt(
        Math.pow(atoms[i].x - atoms[j].x, 2) + Math.pow(atoms[i].y - atoms[j].y, 2)
      );

      if (dist < MIN_DISTANCE && dist > 0) {
        // Move atoms apart slightly
        const dx = atoms[j].x - atoms[i].x;
        const dy = atoms[j].y - atoms[i].y;
        const angle = Math.atan2(dy, dx);

        atoms[i].x -= Math.cos(angle) * (MIN_DISTANCE - dist) * 0.5;
        atoms[i].y -= Math.sin(angle) * (MIN_DISTANCE - dist) * 0.5;
        atoms[j].x += Math.cos(angle) * (MIN_DISTANCE - dist) * 0.5;
        atoms[j].y += Math.sin(angle) * (MIN_DISTANCE - dist) * 0.5;
      }
    }
  }

  return newState;
}

/**
 * Normalize bond orders - ensure bond orders are valid (1, 2, or 3)
 */
function normalizeBondOrders(state: MoleculeState): MoleculeState {
  const newState: MoleculeState = {
    atoms: new Map(state.atoms),
    bonds: new Map(state.bonds),
  };

  newState.bonds.forEach((bond) => {
    if (bond.order < 1) bond.order = 1;
    if (bond.order > 3) bond.order = 3;
    bond.order = Math.round(bond.order);
  });

  return newState;
}

/**
 * Main sanitization function
 */
export function sanitize(state: MoleculeState): MoleculeState {
  let sanitized = { ...state };

  // Apply sanitization steps in order
  sanitized = removeDuplicateBonds(sanitized);
  sanitized = normalizeBondOrders(sanitized);
  sanitized = fixCharges(sanitized);
  sanitized = cleanupGeometry(sanitized);
  sanitized = normalizeHydrogens(sanitized);

  return sanitized;
}
