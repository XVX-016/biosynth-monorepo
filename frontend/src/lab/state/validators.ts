/**
 * State Validators - Runtime invariants and schema checks
 * 
 * Validates molecule state for correctness and can attempt repairs.
 */

import type { MoleculeState } from '../engines/MoleculeStateEngine';

export interface ValidationResult {
  valid: boolean;
  errors: ValidationError[];
  warnings: ValidationWarning[];
  repaired: boolean;
}

export interface ValidationError {
  type: string;
  message: string;
  severity: 'error' | 'warning' | 'info';
  data?: any;
}

export interface ValidationWarning {
  type: string;
  message: string;
  data?: any;
}

/**
 * Validate molecule state for invariants
 */
export function validateState(state: MoleculeState): ValidationResult {
  const errors: ValidationError[] = [];
  const warnings: ValidationWarning[] = [];

  // Check atom invariants
  for (const [atomId, atom] of state.atoms.entries()) {
    // Atom must have valid element
    if (!atom.element || atom.element.trim() === '') {
      errors.push({
        type: 'invalid_atom_element',
        message: `Atom ${atomId} has invalid element`,
        severity: 'error',
        data: { atomId, atom },
      });
    }

    // Atom must have valid position
    if (isNaN(atom.x) || isNaN(atom.y)) {
      errors.push({
        type: 'invalid_atom_position',
        message: `Atom ${atomId} has invalid position`,
        severity: 'error',
        data: { atomId, atom },
      });
    }

    // Charge should be reasonable
    if (Math.abs(atom.charge) > 10) {
      warnings.push({
        type: 'unusual_charge',
        message: `Atom ${atomId} has unusual charge: ${atom.charge}`,
        data: { atomId, charge: atom.charge },
      });
    }
  }

  // Check bond invariants
  for (const [bondId, bond] of state.bonds.entries()) {
    // Bond must reference two different atoms
    if (bond.atoms[0] === bond.atoms[1]) {
      errors.push({
        type: 'self_bond',
        message: `Bond ${bondId} connects atom to itself`,
        severity: 'error',
        data: { bondId, bond },
      });
    }

    // Both atoms must exist
    if (!state.atoms.has(bond.atoms[0])) {
      errors.push({
        type: 'missing_atom',
        message: `Bond ${bondId} references missing atom ${bond.atoms[0]}`,
        severity: 'error',
        data: { bondId, atomId: bond.atoms[0] },
      });
    }

    if (!state.atoms.has(bond.atoms[1])) {
      errors.push({
        type: 'missing_atom',
        message: `Bond ${bondId} references missing atom ${bond.atoms[1]}`,
        severity: 'error',
        data: { bondId, atomId: bond.atoms[1] },
      });
    }

    // Bond order should be valid
    if (bond.order < 1 || bond.order > 3 || !Number.isInteger(bond.order)) {
      errors.push({
        type: 'invalid_bond_order',
        message: `Bond ${bondId} has invalid order: ${bond.order}`,
        severity: 'error',
        data: { bondId, order: bond.order },
      });
    }
  }

  // Check for duplicate bonds
  const bondPairs = new Set<string>();
  for (const [bondId, bond] of state.bonds.entries()) {
    const pair = [bond.atoms[0], bond.atoms[1]].sort().join('_');
    if (bondPairs.has(pair)) {
      warnings.push({
        type: 'duplicate_bond',
        message: `Duplicate bond detected for atoms ${bond.atoms[0]} and ${bond.atoms[1]}`,
        data: { bondId, atoms: bond.atoms },
      });
    }
    bondPairs.add(pair);
  }

  // Check for orphaned atoms (atoms with no bonds)
  for (const [atomId, atom] of state.atoms.entries()) {
    let hasBond = false;
    for (const bond of state.bonds.values()) {
      if (bond.atoms[0] === atomId || bond.atoms[1] === atomId) {
        hasBond = true;
        break;
      }
    }
    if (!hasBond && state.atoms.size > 1) {
      warnings.push({
        type: 'orphaned_atom',
        message: `Atom ${atomId} has no bonds`,
        data: { atomId },
      });
    }
  }

  const valid = errors.filter(e => e.severity === 'error').length === 0;

  return {
    valid,
    errors,
    warnings,
    repaired: false,
  };
}

/**
 * Attempt to repair state issues
 */
export function repairState(state: MoleculeState): { state: MoleculeState; repaired: boolean; fixes: string[] } {
  const fixes: string[] = [];
  let repaired = false;

  // Create a copy to modify
  const repairedState: MoleculeState = {
    atoms: new Map(state.atoms),
    bonds: new Map(state.bonds),
  };

  // Fix missing atoms in bonds
  for (const [bondId, bond] of repairedState.bonds.entries()) {
    if (!repairedState.atoms.has(bond.atoms[0]) || !repairedState.atoms.has(bond.atoms[1])) {
      repairedState.bonds.delete(bondId);
      fixes.push(`Removed bond ${bondId} with missing atoms`);
      repaired = true;
    }
  }

  // Fix self-bonds
  for (const [bondId, bond] of repairedState.bonds.entries()) {
    if (bond.atoms[0] === bond.atoms[1]) {
      repairedState.bonds.delete(bondId);
      fixes.push(`Removed self-bond ${bondId}`);
      repaired = true;
    }
  }

  // Fix invalid bond orders
  for (const [bondId, bond] of repairedState.bonds.entries()) {
    if (bond.order < 1 || bond.order > 3) {
      bond.order = Math.max(1, Math.min(3, Math.round(bond.order)));
      fixes.push(`Fixed bond ${bondId} order to ${bond.order}`);
      repaired = true;
    }
  }

  // Fix invalid atom positions
  for (const [atomId, atom] of repairedState.atoms.entries()) {
    if (isNaN(atom.x) || isNaN(atom.y)) {
      atom.x = atom.x || 0;
      atom.y = atom.y || 0;
      fixes.push(`Fixed atom ${atomId} position`);
      repaired = true;
    }
  }

  // Fix invalid charges
  for (const [atomId, atom] of repairedState.atoms.entries()) {
    if (isNaN(atom.charge)) {
      atom.charge = 0;
      fixes.push(`Fixed atom ${atomId} charge`);
      repaired = true;
    }
  }

  return {
    state: repairedState,
    repaired,
    fixes,
  };
}

/**
 * Validate state and optionally repair
 */
export function validateAndRepair(
  state: MoleculeState,
  autoRepair: boolean = false
): ValidationResult & { repairedState?: MoleculeState; fixes?: string[] } {
  const validation = validateState(state);

  if (!validation.valid && autoRepair) {
    const repair = repairState(state);
    if (repair.repaired) {
      // Re-validate repaired state
      const revalidation = validateState(repair.state);
      return {
        ...revalidation,
        repaired: true,
        repairedState: repair.state,
        fixes: repair.fixes,
      };
    }
  }

  return validation;
}

