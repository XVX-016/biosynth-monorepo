/**
 * AromaticityDetector - Detect aromatic rings in molecules
 */

import type { MoleculeState, Atom, Bond } from '../MoleculeStateEngine';

export interface AromaticRing {
  atoms: string[];
  bonds: string[];
  electronCount: number;
  huckelRule: boolean; // 4n+2 rule
}

export interface AromaticityResult {
  rings: AromaticRing[];
  aromaticAtoms: Set<string>;
  aromaticBonds: Set<string>;
}

/**
 * Detect all rings in molecule
 */
function detectRings(state: MoleculeState): string[][] {
  const rings: string[][] = [];
  const visited = new Set<string>();
  
  const dfs = (
    start: string,
    current: string,
    path: string[],
    visitedInPath: Set<string>
  ) => {
    if (current === start && path.length > 2) {
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
  
  state.atoms.forEach((atom) => {
    if (!visited.has(atom.id)) {
      dfs(atom.id, atom.id, [], new Set());
      visited.add(atom.id);
    }
  });
  
  // Remove duplicates and keep smallest rings
  const uniqueRings = rings.filter((ring, index, self) => {
    const normalized = ring.sort().join('-');
    const firstIndex = self.findIndex((r) => r.sort().join('-') === normalized);
    return index === firstIndex;
  });
  
  // Sort by size and remove supersets
  uniqueRings.sort((a, b) => a.length - b.length);
  const minimalRings: string[][] = [];
  
  uniqueRings.forEach((ring) => {
    const isSubset = minimalRings.some((existing) => {
      const existingSet = new Set(existing);
      return ring.every((atomId) => existingSet.has(atomId));
    });
    
    if (!isSubset) {
      minimalRings.push(ring);
    }
  });
  
  return minimalRings;
}

/**
 * Check if ring is planar (simplified check)
 */
function isPlanar(ring: string[], state: MoleculeState): boolean {
  // For 6-membered rings, assume planar if all atoms are sp2 hybridized
  if (ring.length === 6) {
    return ring.every((atomId) => {
      const atom = state.atoms.get(atomId);
      if (!atom) return false;
      
      const bonds = Array.from(state.bonds.values()).filter(
        (bond) => bond.atoms[0] === atomId || bond.atoms[1] === atomId
      );
      
      const bondOrderSum = bonds.reduce((sum, bond) => sum + bond.order, 0);
      return bondOrderSum >= 3; // sp2 or sp
    });
  }
  
  // For 5-membered rings, also check
  if (ring.length === 5) {
    return true; // Assume planar for common aromatic rings
  }
  
  return false;
}

/**
 * Count pi electrons in ring
 */
function countPiElectrons(ring: string[], state: MoleculeState): number {
  let piElectrons = 0;
  
  ring.forEach((atomId) => {
    const atom = state.atoms.get(atomId);
    if (!atom) return;
    
    const bonds = Array.from(state.bonds.values()).filter(
      (bond) => bond.atoms[0] === atomId || bond.atoms[1] === atomId
    );
    
    // Count double/triple bonds (pi electrons)
    bonds.forEach((bond) => {
      if (bond.order === 2) {
        piElectrons += 2; // Double bond = 2 pi electrons
      } else if (bond.order === 3) {
        piElectrons += 4; // Triple bond = 4 pi electrons (2 pi bonds)
      }
    });
    
    // Lone pairs on heteroatoms contribute
    if (atom.element === 'N' || atom.element === 'O' || atom.element === 'S') {
      const bondOrderSum = bonds.reduce((sum, bond) => sum + bond.order, 0);
      if (bondOrderSum === 2) {
        piElectrons += 2; // Lone pair contributes
      }
    }
  });
  
  return piElectrons;
}

/**
 * Check Huckel's rule (4n+2 pi electrons)
 */
function satisfiesHuckelRule(electronCount: number): boolean {
  // 4n + 2 = electronCount
  // n = (electronCount - 2) / 4
  const n = (electronCount - 2) / 4;
  return n >= 0 && Number.isInteger(n);
}

/**
 * Detect aromaticity in molecule
 */
export function detectAromaticity(state: MoleculeState): AromaticityResult {
  const rings = detectRings(state);
  const aromaticRings: AromaticRing[] = [];
  const aromaticAtoms = new Set<string>();
  const aromaticBonds = new Set<string>();
  
  rings.forEach((ringAtoms) => {
    // Check if ring is planar
    if (!isPlanar(ringAtoms, state)) return;
    
    // Get bonds in ring
    const ringBonds: string[] = [];
    for (let i = 0; i < ringAtoms.length; i++) {
      const atom1 = ringAtoms[i];
      const atom2 = ringAtoms[(i + 1) % ringAtoms.length];
      
      const bond = Array.from(state.bonds.values()).find(
        (b) =>
          (b.atoms[0] === atom1 && b.atoms[1] === atom2) ||
          (b.atoms[0] === atom2 && b.atoms[1] === atom1)
      );
      
      if (bond) {
        ringBonds.push(bond.id);
      }
    }
    
    // Count pi electrons
    const electronCount = countPiElectrons(ringAtoms, state);
    const huckelRule = satisfiesHuckelRule(electronCount);
    
    // Common aromatic rings: benzene (6 electrons), pyridine (6), etc.
    if (huckelRule && electronCount >= 2) {
      aromaticRings.push({
        atoms: ringAtoms,
        bonds: ringBonds,
        electronCount,
        huckelRule: true,
      });
      
      ringAtoms.forEach((atomId) => aromaticAtoms.add(atomId));
      ringBonds.forEach((bondId) => aromaticBonds.add(bondId));
    }
  });
  
  return {
    rings: aromaticRings,
    aromaticAtoms,
    aromaticBonds,
  };
}

/**
 * Main AromaticityDetector class
 */
export class AromaticityDetector {
  /**
   * Detect aromaticity
   */
  detectAromaticity(state: MoleculeState): AromaticityResult {
    return detectAromaticity(state);
  }
  
  /**
   * Check if atom is aromatic
   */
  isAromatic(atomId: string, state: MoleculeState): boolean {
    const result = this.detectAromaticity(state);
    return result.aromaticAtoms.has(atomId);
  }
  
  /**
   * Check if bond is aromatic
   */
  isAromaticBond(bondId: string, state: MoleculeState): boolean {
    const result = this.detectAromaticity(state);
    return result.aromaticBonds.has(bondId);
  }
}

export const aromaticityDetector = new AromaticityDetector();

