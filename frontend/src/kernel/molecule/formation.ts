import { MoleculeGraph } from '@biosynth/engine';
import { ForceField } from '@biosynth/engine';
import { autoBondNewAtom } from '@biosynth/engine';
import { allowedAdditionalBonds, getValence, type ElementSymbol } from '@biosynth/engine';

/**
 * Auto-bond two atoms if they are within bonding distance
 */
export function autoBond(atomA: string, atomB: string, molecule: MoleculeGraph): string | null {
  const atom1 = molecule.atoms.get(atomA);
  const atom2 = molecule.atoms.get(atomB);
  
  if (!atom1 || !atom2) return null;
  
  // Check if already bonded
  const existingBond = Array.from(molecule.bonds.values()).find(
    (b) => (b.a1 === atomA && b.a2 === atomB) || (b.a1 === atomB && b.a2 === atomA)
  );
  if (existingBond) return existingBond.id;
  
  // Check valence constraints
  const bondsA = molecule.getBondsForAtom(atomA);
  const bondsB = molecule.getBondsForAtom(atomB);
  const orderSumA = bondsA.reduce((sum, b) => sum + b.order, 0);
  const orderSumB = bondsB.reduce((sum, b) => sum + b.order, 0);
  
  const remainingA = allowedAdditionalBonds(atom1.element as ElementSymbol, orderSumA);
  const remainingB = allowedAdditionalBonds(atom2.element as ElementSymbol, orderSumB);
  
  if (remainingA <= 0 || remainingB <= 0) return null;
  
  // Calculate distance
  const dx = atom1.position[0] - atom2.position[0];
  const dy = atom1.position[1] - atom2.position[1];
  const dz = atom1.position[2] - atom2.position[2];
  const distance = Math.sqrt(dx * dx + dy * dy + dz * dz);
  
  // Covalent radii (approximate)
  const COVALENT_RADII: Record<string, number> = {
    H: 0.31, C: 0.76, N: 0.71, O: 0.66, F: 0.57,
    P: 1.07, S: 1.05, Cl: 1.02, Br: 1.20, I: 1.39,
  };
  
  const r1 = COVALENT_RADII[atom1.element] ?? 0.8;
  const r2 = COVALENT_RADII[atom2.element] ?? 0.8;
  const maxDist = r1 + r2 + 0.45; // tolerance
  
  if (distance > maxDist) return null;
  
  // Create single bond
  return molecule.addBond(atomA, atomB, 1);
}

/**
 * Validate valency for an atom
 * Returns true if atom has valid valency, false otherwise
 */
export function validateValency(atomId: string, molecule: MoleculeGraph): {
  valid: boolean;
  current: number;
  max: number;
  remaining: number;
} {
  const atom = molecule.atoms.get(atomId);
  if (!atom) {
    return { valid: false, current: 0, max: 0, remaining: 0 };
  }
  
  const bonds = molecule.getBondsForAtom(atomId);
  const currentOrderSum = bonds.reduce((sum, b) => sum + b.order, 0);
  const maxValence = getValence(atom.element as ElementSymbol);
  const remaining = allowedAdditionalBonds(atom.element as ElementSymbol, currentOrderSum);
  const valid = currentOrderSum <= maxValence;
  
  return { valid, current: currentOrderSum, max: maxValence, remaining };
}

/**
 * Recalculate geometry using force field optimization
 */
export function recalcGeometry(
  molecule: MoleculeGraph,
  iterations: number = 20,
  step: number = 0.005
): MoleculeGraph {
  const cloned = molecule.clone();
  ForceField.optimizeGeometry(cloned, iterations, step);
  return cloned;
}

/**
 * Detect rings in the molecule using DFS
 * Returns array of ring atom IDs
 */
export function detectRings(molecule: MoleculeGraph): string[][] {
  const rings: string[][] = [];
  const visited = new Set<string>();
  const path: string[] = [];
  
  function dfs(atomId: string, parentId: string | null) {
    if (visited.has(atomId)) {
      // Found a cycle
      const cycleStart = path.indexOf(atomId);
      if (cycleStart !== -1) {
        const ring = path.slice(cycleStart).concat(atomId);
        // Check if this ring is unique (not a subset of another)
        const isUnique = !rings.some((existingRing) => {
          const existingSet = new Set(existingRing);
          return ring.every((id) => existingSet.has(id));
        });
        if (isUnique && ring.length >= 3) {
          rings.push(ring);
        }
      }
      return;
    }
    
    visited.add(atomId);
    path.push(atomId);
    
    const neighbors = molecule.getNeighbors(atomId);
    for (const neighborId of neighbors) {
      if (neighborId !== parentId) {
        dfs(neighborId, atomId);
      }
    }
    
    path.pop();
    visited.delete(atomId);
  }
  
  // Start DFS from each unvisited atom
  for (const atomId of molecule.atoms.keys()) {
    if (!visited.has(atomId)) {
      dfs(atomId, null);
    }
  }
  
  return rings;
}

/**
 * Merge two molecule fragments
 * Returns new molecule with both fragments combined
 */
export function mergeMoleculeFragments(
  fragment1: MoleculeGraph,
  fragment2: MoleculeGraph,
  offset: [number, number, number] = [0, 0, 0]
): MoleculeGraph {
  const merged = fragment1.clone();
  
  // Map old IDs to new IDs for fragment2
  const idMap = new Map<string, string>();
  
  // Add atoms from fragment2 with offset
  for (const atom of fragment2.atoms.values()) {
    const newPosition: [number, number, number] = [
      atom.position[0] + offset[0],
      atom.position[1] + offset[1],
      atom.position[2] + offset[2],
    ];
    const newId = merged.addAtom({
      element: atom.element,
      position: newPosition,
    });
    idMap.set(atom.id, newId);
  }
  
  // Add bonds from fragment2 using mapped IDs
  for (const bond of fragment2.bonds.values()) {
    const newA1 = idMap.get(bond.a1);
    const newA2 = idMap.get(bond.a2);
    if (newA1 && newA2) {
      merged.addBond(newA1, newA2, bond.order);
    }
  }
  
  return merged;
}

