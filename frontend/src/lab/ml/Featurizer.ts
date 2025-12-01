/**
 * Featurizer - Molecular featurization (ECFP, Graph representation)
 */

import type { MoleculeState, Atom, Bond } from '../engines/MoleculeStateEngine';

export interface GraphFeature {
  nodes: number[][]; // Node features [num_nodes, node_feature_dim]
  edges: number[][]; // Edge indices [num_edges, 2]
  edgeFeatures: number[][]; // Edge features [num_edges, edge_feature_dim]
  nodePositions?: number[][]; // Optional 3D positions [num_nodes, 3]
}

export interface ECFPFeature {
  fingerprint: number[]; // Binary or count fingerprint
  radius: number;
  bitLength: number;
}

/**
 * Element to atomic number mapping
 */
const ELEMENT_TO_ATOMIC_NUMBER: Record<string, number> = {
  H: 1,
  He: 2,
  Li: 3,
  Be: 4,
  B: 5,
  C: 6,
  N: 7,
  O: 8,
  F: 9,
  Ne: 10,
  Na: 11,
  Mg: 12,
  Al: 13,
  Si: 14,
  P: 15,
  S: 16,
  Cl: 17,
  Ar: 18,
  K: 19,
  Ca: 20,
  Fe: 26,
  Cu: 29,
  Zn: 30,
  Br: 35,
  I: 53,
};

/**
 * Compute node features for an atom
 */
function computeNodeFeatures(atom: Atom, bonds: Bond[]): number[] {
  const features: number[] = [];

  // Atomic number (one-hot or embedding index)
  const atomicNum = ELEMENT_TO_ATOMIC_NUMBER[atom.element] || 0;
  features.push(atomicNum);

  // Degree (number of bonds)
  features.push(bonds.length);

  // Bond order sum
  const bondOrderSum = bonds.reduce((sum, bond) => sum + bond.order, 0);
  features.push(bondOrderSum);

  // Formal charge
  features.push(atom.charge || 0);

  // Hybridization estimate (based on bond count and order)
  const hybridization = bondOrderSum <= 2 ? 0 : bondOrderSum <= 3 ? 1 : 2; // sp, sp2, sp3
  features.push(hybridization);

  // Valence electrons (simplified)
  const VALENCE: Record<string, number> = {
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
  };
  features.push(VALENCE[atom.element] || 4);

  // Position (normalized)
  features.push(atom.x / 100); // Normalize
  features.push(atom.y / 100);
  features.push((atom.z || 0) / 100);

  return features;
}

/**
 * Compute edge features for a bond
 */
function computeEdgeFeatures(bond: Bond): number[] {
  const features: number[] = [];

  // Bond order
  features.push(bond.order);

  // Bond type encoding
  features.push(bond.order === 1 ? 1 : 0); // Single
  features.push(bond.order === 2 ? 1 : 0); // Double
  features.push(bond.order === 3 ? 1 : 0); // Triple

  return features;
}

/**
 * Compute graph representation for GNN
 */
export function computeGraph(state: MoleculeState): GraphFeature {
  const atoms = Array.from(state.atoms.values());
  const bonds = Array.from(state.bonds.values());

  // Node features
  const nodes: number[][] = atoms.map((atom) => {
    const atomBonds = bonds.filter(
      (bond) => bond.atoms[0] === atom.id || bond.atoms[1] === atom.id
    );
    return computeNodeFeatures(atom, atomBonds);
  });

  // Edge indices and features
  const edges: number[][] = [];
  const edgeFeatures: number[][] = [];

  bonds.forEach((bond) => {
    const atom1Idx = atoms.findIndex((a) => a.id === bond.atoms[0]);
    const atom2Idx = atoms.findIndex((a) => a.id === bond.atoms[1]);

    if (atom1Idx >= 0 && atom2Idx >= 0) {
      edges.push([atom1Idx, atom2Idx]);
      edgeFeatures.push(computeEdgeFeatures(bond));
    }
  });

  // Node positions (optional, for geometric GNNs)
  const nodePositions = atoms.map((atom) => [
    atom.x / 100,
    atom.y / 100,
    (atom.z || 0) / 100,
  ]);

  return {
    nodes,
    edges,
    edgeFeatures,
    nodePositions,
  };
}

/**
 * Compute ECFP (Extended Connectivity Fingerprint)
 * Simplified version - full ECFP requires RDKit
 */
export function computeECFP(
  state: MoleculeState,
  radius: number = 2,
  bitLength: number = 2048
): ECFPFeature {
  // TODO: Full ECFP implementation requires backend RDKit
  // This is a simplified hash-based fingerprint

  const atoms = Array.from(state.atoms.values());
  const bonds = Array.from(state.bonds.values());
  const fingerprint = new Array(bitLength).fill(0);

  // Simple hash-based fingerprint
  atoms.forEach((atom) => {
    // Hash atom type
    const atomHash = hashString(atom.element) % bitLength;
    fingerprint[atomHash] = 1;

    // Hash atom + neighbors (radius 1)
    const neighbors = bonds
      .filter((bond) => bond.atoms[0] === atom.id || bond.atoms[1] === atom.id)
      .map((bond) => {
        const otherId = bond.atoms[0] === atom.id ? bond.atoms[1] : bond.atoms[0];
        return state.atoms.get(otherId);
      })
      .filter((a): a is Atom => a !== undefined);

    neighbors.forEach((neighbor) => {
      const pairHash = hashString(`${atom.element}-${neighbor.element}`) % bitLength;
      fingerprint[pairHash] = 1;
    });
  });

  return {
    fingerprint,
    radius,
    bitLength,
  };
}

/**
 * Simple string hash function
 */
function hashString(str: string): number {
  let hash = 0;
  for (let i = 0; i < str.length; i++) {
    const char = str.charCodeAt(i);
    hash = (hash << 5) - hash + char;
    hash = hash & hash; // Convert to 32-bit integer
  }
  return Math.abs(hash);
}

/**
 * Main Featurizer class
 */
export class Featurizer {
  /**
   * Compute graph features for GNN
   */
  computeGraph(state: MoleculeState): GraphFeature {
    return computeGraph(state);
  }

  /**
   * Compute ECFP fingerprint
   */
  computeECFP(state: MoleculeState, radius?: number, bitLength?: number): ECFPFeature {
    return computeECFP(state, radius, bitLength);
  }

  /**
   * Compute both graph and ECFP features
   */
  computeAll(state: MoleculeState): {
    graph: GraphFeature;
    ecfp: ECFPFeature;
  } {
    return {
      graph: this.computeGraph(state),
      ecfp: this.computeECFP(state),
    };
  }
}

export const featurizer = new Featurizer();

