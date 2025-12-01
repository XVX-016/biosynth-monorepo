/**
 * Attention mapping utilities
 * 
 * Maps backend attention weights to frontend bond/atom IDs
 */

export interface AttentionMapping {
  bondId?: string;
  uIdx: number;
  vIdx: number;
  uAtomId: string;
  vAtomId: string;
  attention: number;
}

/**
 * Map edge_index to bond IDs
 */
export function mapAttentionToBonds(
  edgeIndex: number[][],
  feAtomOrder: string[],
  bondsByAtomPair: Record<string, string>,
  attentions: number[]
): Map<string, number> {
  /**
   * edgeIndex: [2][E] - backend edge indices
   * feAtomOrder: [atomId0, atomId1, ...] - frontend atom IDs in order
   * bondsByAtomPair: {"a0_a1": "bond123", ...} - bond ID lookup
   * attentions: [E] - attention weights per edge
   */
  
  const bondAttentions = new Map<string, number>();
  const E = edgeIndex[0]?.length || 0;
  
  for (let e = 0; e < E && e < attentions.length; e++) {
    const u = edgeIndex[0][e];
    const v = edgeIndex[1][e];
    
    if (u >= feAtomOrder.length || v >= feAtomOrder.length) {
      continue;
    }
    
    const uAtomId = feAtomOrder[u];
    const vAtomId = feAtomOrder[v];
    
    if (!uAtomId || !vAtomId) {
      continue;
    }
    
    // Try both orderings
    const key1 = `${uAtomId}_${vAtomId}`;
    const key2 = `${vAtomId}_${uAtomId}`;
    const bondId = bondsByAtomPair[key1] || bondsByAtomPair[key2];
    
    if (bondId) {
      bondAttentions.set(bondId, attentions[e]);
    }
  }
  
  return bondAttentions;
}

/**
 * Map node importance to atom IDs
 */
export function mapNodeImportanceToAtoms(
  nodeImportance: number[],
  feAtomOrder: string[]
): Map<string, number> {
  const atomImportance = new Map<string, number>();
  
  nodeImportance.forEach((importance, idx) => {
    if (idx < feAtomOrder.length) {
      const atomId = feAtomOrder[idx];
      if (atomId) {
        atomImportance.set(atomId, importance);
      }
    }
  });
  
  return atomImportance;
}

/**
 * Create bondsByAtomPair lookup from molecule state
 */
export function createBondLookup(
  atoms: Array<{ id: string }>,
  bonds: Array<{ id: string; atoms: [string, string] }>
): Record<string, string> {
  const lookup: Record<string, string> = {};
  
  bonds.forEach((bond) => {
    const [a1, a2] = bond.atoms;
    // Store both orderings
    lookup[`${a1}_${a2}`] = bond.id;
    lookup[`${a2}_${a1}`] = bond.id;
  });
  
  return lookup;
}

/**
 * Get atom order from molecule state
 */
export function getAtomOrder(
  atoms: Array<{ id: string }>
): string[] {
  return atoms.map((atom) => atom.id);
}

/**
 * Complete attention mapping workflow
 */
export function mapAttentionWeights(
  edgeIndex: number[][],
  attentions: number[],
  nodeImportance: number[] | undefined,
  atoms: Array<{ id: string }>,
  bonds: Array<{ id: string; atoms: [string, string] }>
): {
  bondAttentions: Map<string, number>;
  atomImportance: Map<string, number>;
} {
  const atomOrder = getAtomOrder(atoms);
  const bondLookup = createBondLookup(atoms, bonds);
  
  const bondAttentions = mapAttentionToBonds(
    edgeIndex,
    atomOrder,
    bondLookup,
    attentions
  );
  
  const atomImportanceMap = nodeImportance
    ? mapNodeImportanceToAtoms(nodeImportance, atomOrder)
    : new Map<string, number>();
  
  return {
    bondAttentions,
    atomImportance: atomImportanceMap,
  };
}

