/**
 * AttentionVisualizer - Visualize attention weights on molecule
 * 
 * This component applies attention weights to bonds/atoms in the Editor2D
 */

import React, { useEffect } from 'react';
import { useLab } from '../hooks/useLab';
import type { MoleculeState } from '../engines/MoleculeStateEngine';

export interface AttentionVisualizerProps {
  attentions: {
    edgeAttentions: number[];
    edgeIndex: number[][];
    nodeImportance?: number[];
  };
  molecule: MoleculeState;
  layerIndex?: number;
}

/**
 * Map attention weights to bond IDs
 */
function mapAttentionsToBonds(
  attentions: AttentionVisualizerProps['attentions'],
  molecule: MoleculeState
): Map<string, number> {
  const bondAttentions = new Map<string, number>();
  const atoms = Array.from(molecule.atoms.values());

  // Create atom index to atom ID mapping
  const atomIndexToId = new Map<number, string>();
  atoms.forEach((atom, idx) => {
    atomIndexToId.set(idx, atom.id);
  });

  // Map edge indices to bond IDs
  attentions.edgeIndex.forEach((edge, idx) => {
    const [u, v] = edge;
    const atom1Id = atomIndexToId.get(u);
    const atom2Id = atomIndexToId.get(v);

    if (atom1Id && atom2Id) {
      // Find bond connecting these atoms
      const bond = Array.from(molecule.bonds.values()).find(
        (b) =>
          (b.atoms[0] === atom1Id && b.atoms[1] === atom2Id) ||
          (b.atoms[0] === atom2Id && b.atoms[1] === atom1Id)
      );

      if (bond && idx < attentions.edgeAttentions.length) {
        bondAttentions.set(bond.id, attentions.edgeAttentions[idx]);
      }
    }
  });

  return bondAttentions;
}

/**
 * Map node importance to atom IDs
 */
function mapNodeImportanceToAtoms(
  nodeImportance: number[],
  molecule: MoleculeState
): Map<string, number> {
  const atomImportance = new Map<string, number>();
  const atoms = Array.from(molecule.atoms.values());

  nodeImportance.forEach((importance, idx) => {
    if (idx < atoms.length) {
      atomImportance.set(atoms[idx].id, importance);
    }
  });

  return atomImportance;
}

/**
 * AttentionVisualizer component
 * 
 * This component doesn't render UI directly - it provides data
 * for Editor2D to visualize attention weights
 */
export default function AttentionVisualizer({
  attentions,
  molecule,
}: AttentionVisualizerProps) {
  // This component is used to compute attention mappings
  // The actual visualization happens in Editor2D
  
  const bondAttentions = mapAttentionsToBonds(attentions, molecule);
  const atomImportance = attentions.nodeImportance
    ? mapNodeImportanceToAtoms(attentions.nodeImportance, molecule)
    : new Map<string, number>();

  // Export via custom hook or context
  // For now, this is a utility component
  
  return null; // No direct rendering
}

/**
 * Hook to get attention weights for visualization
 */
export function useAttentionWeights(
  attentions: AttentionVisualizerProps['attentions'] | null,
  molecule: MoleculeState | null
): {
  bondAttentions: Map<string, number>;
  atomImportance: Map<string, number>;
} {
  if (!attentions || !molecule) {
    return {
      bondAttentions: new Map(),
      atomImportance: new Map(),
    };
  }

  const bondAttentions = mapAttentionsToBonds(attentions, molecule);
  const atomImportance = attentions.nodeImportance
    ? mapNodeImportanceToAtoms(attentions.nodeImportance, molecule)
    : new Map<string, number>();

  return { bondAttentions, atomImportance };
}

