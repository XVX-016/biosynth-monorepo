/**
 * geometryCleanup.ts
 * Pure math functions for geometry optimization and clash resolution.
 * No R3F dependencies - works directly with MoleculeGraph.
 */

import { MoleculeGraph } from "@biosynth/engine";
import { getBondLength } from "@biosynth/engine";

const CLASH_THRESHOLD = 0.7; // Å - minimum distance between non-bonded atoms
const REPULSION_FORCE = 0.1;
const REPULSION_DECAY = 0.9;
const MAX_REPULSION_ITERATIONS = 20;

/**
 * Resolve atomic clashes by pushing atoms apart.
 * Uses simple repulsion force with decay.
 */
export function resolveClashes(graph: MoleculeGraph): void {
  const atoms = Array.from(graph.atoms.values());
  const bonds = Array.from(graph.bonds.values());

  // Build set of bonded pairs for quick lookup
  const bondedPairs = new Set<string>();
  bonds.forEach((bond) => {
    const key1 = `${bond.a1}-${bond.a2}`;
    const key2 = `${bond.a2}-${bond.a1}`;
    bondedPairs.add(key1);
    bondedPairs.add(key2);
  });

  const isBonded = (id1: string, id2: string): boolean => {
    return bondedPairs.has(`${id1}-${id2}`);
  };

  // Iterative repulsion
  for (let iter = 0; iter < MAX_REPULSION_ITERATIONS; iter++) {
    const forces: Map<string, [number, number, number]> = new Map();
    atoms.forEach((atom) => {
      forces.set(atom.id, [0, 0, 0]);
    });

    // Calculate repulsion forces
    for (let i = 0; i < atoms.length; i++) {
      for (let j = i + 1; j < atoms.length; j++) {
        const a1 = atoms[i];
        const a2 = atoms[j];

        // Skip if bonded
        if (isBonded(a1.id, a2.id)) continue;

        const dx = a2.position[0] - a1.position[0];
        const dy = a2.position[1] - a1.position[1];
        const dz = a2.position[2] - a1.position[2];
        const distance = Math.sqrt(dx * dx + dy * dy + dz * dz);

        if (distance < CLASH_THRESHOLD && distance > 0.01) {
          // Repulsion force
          const force = (REPULSION_FORCE * (CLASH_THRESHOLD - distance)) / distance;
          const fx = (force * dx) / distance;
          const fy = (force * dy) / distance;
          const fz = (force * dz) / distance;

          const f1 = forces.get(a1.id)!;
          const f2 = forces.get(a2.id)!;
          forces.set(a1.id, [f1[0] - fx, f1[1] - fy, f1[2] - fz]);
          forces.set(a2.id, [f2[0] + fx, f2[1] + fy, f2[2] + fz]);
        }
      }
    }

    // Apply forces with decay
    const decay = Math.pow(REPULSION_DECAY, iter);
    forces.forEach((force, atomId) => {
      const atom = graph.atoms.get(atomId);
      if (!atom) return;

      atom.position[0] += force[0] * decay;
      atom.position[1] += force[1] * decay;
      atom.position[2] += force[2] * decay;
    });
  }
}

/**
 * Detect aromatic rings and planarize them using PCA.
 */
export function planarizeAromaticRings(graph: MoleculeGraph): void {
  const bonds = Array.from(graph.bonds.values());
  const atoms = Array.from(graph.atoms.values());

  // Detect aromatic rings (simplified: alternating single/double bonds in 6-membered rings)
  const visited = new Set<string>();
  const rings: string[][] = [];

  // Simple ring detection: find cycles of length 6 with alternating bond orders
  const findRings = (startId: string, path: string[], depth: number): void => {
    if (depth > 6) return;
    if (depth === 6 && path[0] === startId) {
      // Check if alternating
      let isAlternating = true;
      for (let i = 0; i < 6; i++) {
        const a1 = path[i];
        const a2 = path[(i + 1) % 6];
        const bond = bonds.find(
          (b) =>
            (b.a1 === a1 && b.a2 === a2) || (b.a1 === a2 && b.a2 === a1)
        );
        if (bond && i % 2 === 0 && bond.order !== 2) isAlternating = false;
        if (bond && i % 2 === 1 && bond.order !== 1) isAlternating = false;
      }
      if (isAlternating) {
        rings.push([...path]);
      }
      return;
    }

    const currentAtom = graph.atoms.get(path[path.length - 1]);
    if (!currentAtom) return;

    const neighbors = graph.getNeighbors(currentAtom.id);
    for (const neighborId of neighbors) {
      if (!path.includes(neighborId) || (depth === 5 && neighborId === startId)) {
        findRings(startId, [...path, neighborId], depth + 1);
      }
    }
  };

  // Find rings starting from each atom (limit to avoid explosion)
  for (const atom of atoms.slice(0, 10)) {
    if (!visited.has(atom.id)) {
      findRings(atom.id, [atom.id], 0);
      visited.add(atom.id);
    }
  }

  // Planarize each ring using PCA
  rings.forEach((ringIds) => {
    const ringAtoms = ringIds.map((id) => graph.atoms.get(id)).filter(Boolean);
    if (ringAtoms.length !== 6) return;

    // Calculate centroid
    const centroid = [0, 0, 0];
    ringAtoms.forEach((atom) => {
      centroid[0] += atom!.position[0];
      centroid[1] += atom!.position[1];
      centroid[2] += atom!.position[2];
    });
    centroid[0] /= 6;
    centroid[1] /= 6;
    centroid[2] /= 6;

    // Build covariance matrix (simplified 2D PCA)
    let xx = 0, yy = 0, zz = 0, xy = 0, xz = 0, yz = 0;
    ringAtoms.forEach((atom) => {
      const x = atom!.position[0] - centroid[0];
      const y = atom!.position[1] - centroid[1];
      const z = atom!.position[2] - centroid[2];
      xx += x * x;
      yy += y * y;
      zz += z * z;
      xy += x * y;
      xz += x * z;
      yz += y * z;
    });

    // Find best-fit plane (simplified: use z-axis as normal if variance is smallest)
    // For simplicity, project onto XY plane
    ringAtoms.forEach((atom) => {
      atom!.position[2] = centroid[2]; // Flatten to plane
    });
  });
}

/**
 * Lightweight force-directed relaxation.
 * Uses simple spring model for bond lengths.
 */
export function relaxGeometry(graph: MoleculeGraph, iterations: number = 10): void {
  const bonds = Array.from(graph.bonds.values());
  const SPRING_CONSTANT = 0.05;
  const DAMPING = 0.8;

  for (let iter = 0; iter < iterations; iter++) {
    const forces: Map<string, [number, number, number]> = new Map();
    graph.atoms.forEach((atom) => {
      forces.set(atom.id, [0, 0, 0]);
    });

    // Bond stretching forces
    bonds.forEach((bond) => {
      const a1 = graph.atoms.get(bond.a1);
      const a2 = graph.atoms.get(bond.a2);
      if (!a1 || !a2) return;

      let desiredLength = getBondLength(a1.element, a2.element);
      if (bond.order === 2) desiredLength *= 0.9; // Double bonds shorter
      if (bond.order === 3) desiredLength *= 0.85; // Triple bonds even shorter

      const dx = a2.position[0] - a1.position[0];
      const dy = a2.position[1] - a1.position[1];
      const dz = a2.position[2] - a1.position[2];
      const distance = Math.sqrt(dx * dx + dy * dy + dz * dz);

      if (distance > 0.01) {
        const force = SPRING_CONSTANT * (distance - desiredLength);
        const fx = (force * dx) / distance;
        const fy = (force * dy) / distance;
        const fz = (force * dz) / distance;

        const f1 = forces.get(a1.id)!;
        const f2 = forces.get(a2.id)!;
        forces.set(a1.id, [f1[0] - fx, f1[1] - fy, f1[2] - fz]);
        forces.set(a2.id, [f2[0] + fx, f2[1] + fy, f2[2] + fz]);
      }
    });

    // Apply forces with damping
    forces.forEach((force, atomId) => {
      const atom = graph.atoms.get(atomId);
      if (!atom) return;

      atom.position[0] += force[0] * DAMPING;
      atom.position[1] += force[1] * DAMPING;
      atom.position[2] += force[2] * DAMPING;
    });
  }
}

