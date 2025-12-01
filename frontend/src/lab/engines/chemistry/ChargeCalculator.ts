/**
 * ChargeCalculator - Calculate formal charges for atoms
 */

import type { MoleculeState, Atom, Bond } from '../MoleculeStateEngine';

export interface ChargeResult {
  charges: Map<string, number>;
  totalCharge: number;
}

/**
 * Calculate formal charge for an atom
 * Formal charge = valence electrons - (bonds + lone pairs)
 * Simplified: FC = group number - bonds - lone pairs
 */
function calculateFormalCharge(
  atom: Atom,
  bonds: Bond[],
  lonePairs: number = 0
): number {
  const VALENCE_ELECTRONS: Record<string, number> = {
    H: 1,
    C: 4,
    N: 5,
    O: 6,
    F: 7,
    Cl: 7,
    Br: 7,
    I: 7,
    S: 6,
    P: 5,
    B: 3,
    Si: 4,
  };
  
  const valenceElectrons = VALENCE_ELECTRONS[atom.element] || 4;
  const bondCount = bonds.length;
  const bondOrderSum = bonds.reduce((sum, bond) => sum + bond.order, 0);
  
  // Formal charge = valence electrons - (bonding electrons / 2) - lone pairs
  // Bonding electrons = sum of bond orders
  const bondingElectrons = bondOrderSum;
  const formalCharge = valenceElectrons - bondingElectrons - lonePairs * 2;
  
  return formalCharge;
}

/**
 * Estimate lone pairs based on bonding
 */
function estimateLonePairs(atom: Atom, bonds: Bond[]): number {
  const MAX_BONDS: Record<string, number> = {
    H: 1,
    C: 4,
    N: 4,
    O: 2,
    F: 1,
    Cl: 1,
    Br: 1,
    I: 1,
    S: 6,
    P: 5,
    B: 3,
    Si: 4,
  };
  
  const maxBonds = MAX_BONDS[atom.element] || 4;
  const bondOrderSum = bonds.reduce((sum, bond) => sum + bond.order, 0);
  
  if (bondOrderSum >= maxBonds) {
    return 0; // No lone pairs
  }
  
  // Estimate based on octet rule
  const VALENCE_ELECTRONS: Record<string, number> = {
    H: 2, // Exception: H wants 2 electrons
    C: 8,
    N: 8,
    O: 8,
    F: 8,
    Cl: 8,
    Br: 8,
    I: 8,
    S: 8,
    P: 8,
    B: 6,
    Si: 8,
  };
  
  const targetElectrons = VALENCE_ELECTRONS[atom.element] || 8;
  const currentElectrons = bondOrderSum * 2; // Each bond contributes 2 electrons
  const missingElectrons = targetElectrons - currentElectrons;
  
  return Math.max(0, Math.floor(missingElectrons / 2));
}

/**
 * Compute formal charges for all atoms
 */
export function computeCharges(state: MoleculeState): ChargeResult {
  const charges = new Map<string, number>();
  
  state.atoms.forEach((atom) => {
    const bonds = Array.from(state.bonds.values()).filter(
      (bond) => bond.atoms[0] === atom.id || bond.atoms[1] === atom.id
    );
    
    const lonePairs = estimateLonePairs(atom, bonds);
    const formalCharge = calculateFormalCharge(atom, bonds, lonePairs);
    
    charges.set(atom.id, formalCharge);
  });
  
  const totalCharge = Array.from(charges.values()).reduce((sum, charge) => sum + charge, 0);
  
  return { charges, totalCharge };
}

/**
 * Apply computed charges to molecule state
 */
export function applyCharges(state: MoleculeState, charges: Map<string, number>): void {
  charges.forEach((charge, atomId) => {
    const atom = state.atoms.get(atomId);
    if (atom) {
      atom.charge = charge;
    }
  });
}

/**
 * Main ChargeCalculator class
 */
export class ChargeCalculator {
  /**
   * Compute charges for molecule
   */
  computeCharges(state: MoleculeState): ChargeResult {
    return computeCharges(state);
  }
  
  /**
   * Compute and apply charges
   */
  computeAndApply(state: MoleculeState): ChargeResult {
    const result = this.computeCharges(state);
    applyCharges(state, result.charges);
    return result;
  }
}

export const chargeCalculator = new ChargeCalculator();

