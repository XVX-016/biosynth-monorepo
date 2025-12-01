/**
 * AngleFixer - Fix bond angles and hybridization
 */

import type { MoleculeState, Atom, Bond } from '../MoleculeStateEngine';

export interface AngleFixResult {
  fixed: boolean;
  changes: Array<{ atomId: string; oldAngle: number; newAngle: number }>;
}

/**
 * Calculate ideal bond angle based on hybridization
 */
function getIdealAngle(hybridization: 'sp' | 'sp2' | 'sp3'): number {
  switch (hybridization) {
    case 'sp':
      return Math.PI; // 180 degrees (linear)
    case 'sp2':
      return (2 * Math.PI) / 3; // 120 degrees (trigonal)
    case 'sp3':
      return Math.acos(-1 / 3); // ~109.5 degrees (tetrahedral)
    default:
      return Math.acos(-1 / 3); // Default to tetrahedral
  }
}

/**
 * Determine hybridization from bond count and geometry
 */
function determineHybridization(atom: Atom, bonds: Bond[]): 'sp' | 'sp2' | 'sp3' {
  const bondCount = bonds.length;
  const bondOrderSum = bonds.reduce((sum, bond) => sum + bond.order, 0);
  
  if (bondCount === 2 && bondOrderSum === 2) {
    return 'sp'; // Linear (e.g., CO2)
  } else if (bondCount === 3 && bondOrderSum >= 3) {
    return 'sp2'; // Trigonal (e.g., C=C)
  } else {
    return 'sp3'; // Tetrahedral (default)
  }
}

/**
 * Calculate current bond angle between three atoms
 */
function calculateAngle(
  atom1: Atom,
  center: Atom,
  atom2: Atom
): number {
  const v1 = { x: atom1.x - center.x, y: atom1.y - center.y };
  const v2 = { x: atom2.x - center.x, y: atom2.y - center.y };
  
  const dot = v1.x * v2.x + v1.y * v2.y;
  const mag1 = Math.sqrt(v1.x * v1.x + v1.y * v1.y);
  const mag2 = Math.sqrt(v2.x * v2.x + v2.y * v2.y);
  
  if (mag1 === 0 || mag2 === 0) return Math.PI;
  
  return Math.acos(Math.max(-1, Math.min(1, dot / (mag1 * mag2))));
}

/**
 * Fix bond angles for an atom
 */
function fixAtomAngles(
  atom: Atom,
  bonds: Bond[],
  state: MoleculeState,
  threshold: number = 0.1
): boolean {
  if (bonds.length < 2) return false;
  
  const hybridization = determineHybridization(atom, bonds);
  const idealAngle = getIdealAngle(hybridization);
  
  let fixed = false;
  
  // For each pair of bonds, check and fix angle
  for (let i = 0; i < bonds.length; i++) {
    for (let j = i + 1; j < bonds.length; j++) {
      const bond1 = bonds[i];
      const bond2 = bonds[j];
      
      const otherAtom1Id = bond1.atoms[0] === atom.id ? bond1.atoms[1] : bond1.atoms[0];
      const otherAtom2Id = bond2.atoms[0] === atom.id ? bond2.atoms[1] : bond2.atoms[0];
      
      const otherAtom1 = state.atoms.get(otherAtom1Id);
      const otherAtom2 = state.atoms.get(otherAtom2Id);
      
      if (!otherAtom1 || !otherAtom2) continue;
      
      const currentAngle = calculateAngle(otherAtom1, atom, otherAtom2);
      const angleDiff = Math.abs(currentAngle - idealAngle);
      
      if (angleDiff > threshold) {
        // Adjust positions to fix angle
        const angleToFix = idealAngle - currentAngle;
        
        // Rotate otherAtom2 around atom to fix angle
        const dx = otherAtom2.x - atom.x;
        const dy = otherAtom2.y - atom.y;
        const dist = Math.sqrt(dx * dx + dy * dy);
        
        const currentAngleToAtom1 = Math.atan2(
          otherAtom1.y - atom.y,
          otherAtom1.x - atom.x
        );
        const newAngle = currentAngleToAtom1 + idealAngle;
        
        otherAtom2.x = atom.x + Math.cos(newAngle) * dist;
        otherAtom2.y = atom.y + Math.sin(newAngle) * dist;
        
        fixed = true;
      }
    }
  }
  
  return fixed;
}

/**
 * Fix hybridization by adjusting bond orders
 */
function fixHybridization(
  atom: Atom,
  bonds: Bond[],
  state: MoleculeState
): boolean {
  const hybridization = determineHybridization(atom, bonds);
  const bondCount = bonds.length;
  const bondOrderSum = bonds.reduce((sum, bond) => sum + bond.order, 0);
  
  let fixed = false;
  
  // For carbon: ensure proper hybridization
  if (atom.element === 'C') {
    if (hybridization === 'sp2' && bondCount === 2 && bondOrderSum < 3) {
      // Should have a double bond
      const singleBond = bonds.find((b) => b.order === 1);
      if (singleBond) {
        singleBond.order = 2;
        fixed = true;
      }
    } else if (hybridization === 'sp' && bondCount === 2 && bondOrderSum < 4) {
      // Should have triple bond or two double bonds
      if (bondOrderSum === 2) {
        bonds[0].order = 2;
        bonds[1].order = 2;
        fixed = true;
      }
    }
  }
  
  return fixed;
}

/**
 * Fix all bond angles in molecule
 */
export function fixAngles(
  state: MoleculeState,
  threshold: number = 0.1
): AngleFixResult {
  const changes: AngleFixResult['changes'] = [];
  let fixed = false;
  
  state.atoms.forEach((atom) => {
    const bonds = Array.from(state.bonds.values()).filter(
      (bond) => bond.atoms[0] === atom.id || bond.atoms[1] === atom.id
    );
    
    if (bonds.length >= 2) {
      const oldAngles: number[] = [];
      
      // Calculate old angles
      for (let i = 0; i < bonds.length; i++) {
        for (let j = i + 1; j < bonds.length; j++) {
          const bond1 = bonds[i];
          const bond2 = bonds[j];
          const otherAtom1Id = bond1.atoms[0] === atom.id ? bond1.atoms[1] : bond1.atoms[0];
          const otherAtom2Id = bond2.atoms[0] === atom.id ? bond2.atoms[1] : bond2.atoms[0];
          const otherAtom1 = state.atoms.get(otherAtom1Id);
          const otherAtom2 = state.atoms.get(otherAtom2Id);
          
          if (otherAtom1 && otherAtom2) {
            oldAngles.push(calculateAngle(otherAtom1, atom, otherAtom2));
          }
        }
      }
      
      if (fixAtomAngles(atom, bonds, state, threshold)) {
        fixed = true;
        
        // Calculate new angles
        const newAngles: number[] = [];
        for (let i = 0; i < bonds.length; i++) {
          for (let j = i + 1; j < bonds.length; j++) {
            const bond1 = bonds[i];
            const bond2 = bonds[j];
            const otherAtom1Id = bond1.atoms[0] === atom.id ? bond1.atoms[1] : bond1.atoms[0];
            const otherAtom2Id = bond2.atoms[0] === atom.id ? bond2.atoms[1] : bond2.atoms[0];
            const otherAtom1 = state.atoms.get(otherAtom1Id);
            const otherAtom2 = state.atoms.get(otherAtom2Id);
            
            if (otherAtom1 && otherAtom2) {
              newAngles.push(calculateAngle(otherAtom1, atom, otherAtom2));
            }
          }
        }
        
        if (oldAngles.length > 0 && newAngles.length > 0) {
          changes.push({
            atomId: atom.id,
            oldAngle: oldAngles[0],
            newAngle: newAngles[0],
          });
        }
      }
    }
  });
  
  return { fixed, changes };
}

/**
 * Fix hybridization for all atoms
 */
export function fixHybridization(state: MoleculeState): boolean {
  let fixed = false;
  
  state.atoms.forEach((atom) => {
    const bonds = Array.from(state.bonds.values()).filter(
      (bond) => bond.atoms[0] === atom.id || bond.atoms[1] === atom.id
    );
    
    if (fixHybridization(atom, bonds, state)) {
      fixed = true;
    }
  });
  
  return fixed;
}

/**
 * Main AngleFixer class
 */
export class AngleFixer {
  /**
   * Fix all bond angles
   */
  fixAngles(state: MoleculeState, threshold?: number): AngleFixResult {
    return fixAngles(state, threshold);
  }
  
  /**
   * Fix hybridization
   */
  fixHybridization(state: MoleculeState): boolean {
    return fixHybridization(state);
  }
  
  /**
   * Fix both angles and hybridization
   */
  fixAll(state: MoleculeState, threshold?: number): { angles: AngleFixResult; hybridization: boolean } {
    return {
      angles: this.fixAngles(state, threshold),
      hybridization: this.fixHybridization(state),
    };
  }
}

export const angleFixer = new AngleFixer();

