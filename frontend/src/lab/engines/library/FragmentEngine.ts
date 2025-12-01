/**
 * FragmentEngine - Manage molecular fragments library
 */

import type { MoleculeState, Atom, Bond } from '../MoleculeStateEngine';

export interface Fragment {
  id: string;
  name: string;
  state: MoleculeState;
  attachmentPoints: string[]; // Atom IDs where fragment can attach
  tags: string[];
}

export interface FragmentMatch {
  fragment: Fragment;
  atoms: string[];
  bonds: string[];
  score: number;
}

/**
 * Fragment library storage
 */
class FragmentLibrary {
  private fragments = new Map<string, Fragment>();
  
  /**
   * Add fragment to library
   */
  add(fragment: Fragment): void {
    this.fragments.set(fragment.id, fragment);
  }
  
  /**
   * Get fragment by ID
   */
  get(id: string): Fragment | undefined {
    return this.fragments.get(id);
  }
  
  /**
   * Get all fragments
   */
  getAll(): Fragment[] {
    return Array.from(this.fragments.values());
  }
  
  /**
   * Search fragments by tag
   */
  findByTag(tag: string): Fragment[] {
    return Array.from(this.fragments.values()).filter((f) => f.tags.includes(tag));
  }
  
  /**
   * Search fragments by name
   */
  findByName(name: string): Fragment[] {
    const lowerName = name.toLowerCase();
    return Array.from(this.fragments.values()).filter((f) =>
      f.name.toLowerCase().includes(lowerName)
    );
  }
  
  /**
   * Remove fragment
   */
  remove(id: string): boolean {
    return this.fragments.delete(id);
  }
  
  /**
   * Clear all fragments
   */
  clear(): void {
    this.fragments.clear();
  }
}

/**
 * Match fragment in molecule
 */
function matchFragment(
  state: MoleculeState,
  fragment: Fragment
): FragmentMatch | null {
  // Simple matching: check if fragment atoms exist in molecule
  // Full implementation would use subgraph isomorphism
  
  const fragmentAtoms = Array.from(fragment.state.atoms.values());
  const fragmentBonds = Array.from(fragment.state.bonds.values());
  
  // Try to find matching atoms
  const matchedAtoms: string[] = [];
  const matchedBonds: string[] = [];
  
  // For each fragment atom, try to find matching atom in molecule
  fragmentAtoms.forEach((fragAtom) => {
    const matchingAtom = Array.from(state.atoms.values()).find(
      (atom) => atom.element === fragAtom.element
    );
    
    if (matchingAtom && !matchedAtoms.includes(matchingAtom.id)) {
      matchedAtoms.push(matchingAtom.id);
    }
  });
  
  // Check if we found enough atoms
  if (matchedAtoms.length < fragmentAtoms.length * 0.8) {
    return null; // Not enough matches
  }
  
  // Try to match bonds
  fragmentBonds.forEach((fragBond) => {
    const atom1Idx = fragmentAtoms.findIndex((a) => a.id === fragBond.atoms[0]);
    const atom2Idx = fragmentAtoms.findIndex((a) => a.id === fragBond.atoms[1]);
    
    if (atom1Idx >= 0 && atom2Idx >= 0 && atom1Idx < matchedAtoms.length && atom2Idx < matchedAtoms.length) {
      const molAtom1 = matchedAtoms[atom1Idx];
      const molAtom2 = matchedAtoms[atom2Idx];
      
      const matchingBond = Array.from(state.bonds.values()).find(
        (bond) =>
          ((bond.atoms[0] === molAtom1 && bond.atoms[1] === molAtom2) ||
            (bond.atoms[0] === molAtom2 && bond.atoms[1] === molAtom1)) &&
          bond.order === fragBond.order
      );
      
      if (matchingBond) {
        matchedBonds.push(matchingBond.id);
      }
    }
  });
  
  const score = (matchedAtoms.length / fragmentAtoms.length) * 0.7 +
                (matchedBonds.length / fragmentBonds.length) * 0.3;
  
  if (score < 0.5) {
    return null;
  }
  
  return {
    fragment,
    atoms: matchedAtoms,
    bonds: matchedBonds,
    score,
  };
}

/**
 * Main FragmentEngine class
 */
export class FragmentEngine {
  private library = new FragmentLibrary();
  
  /**
   * Add fragment to library
   */
  addFragment(fragment: Fragment): void {
    this.library.add(fragment);
  }
  
  /**
   * Find fragment by pattern (searches library)
   */
  findFragment(pattern: string | MoleculeState): Fragment[] {
    if (typeof pattern === 'string') {
      // Search by name or tag
      return [
        ...this.library.findByName(pattern),
        ...this.library.findByTag(pattern),
      ];
    } else {
      // Search by matching structure
      const matches: FragmentMatch[] = [];
      
      this.library.getAll().forEach((fragment) => {
        const match = matchFragment(pattern, fragment);
        if (match) {
          matches.push(match);
        }
      });
      
      // Sort by score
      matches.sort((a, b) => b.score - a.score);
      
      return matches.map((m) => m.fragment);
    }
  }
  
  /**
   * Get fragment by ID
   */
  getFragment(id: string): Fragment | undefined {
    return this.library.get(id);
  }
  
  /**
   * Get all fragments
   */
  getAllFragments(): Fragment[] {
    return this.library.getAll();
  }
  
  /**
   * Remove fragment
   */
  removeFragment(id: string): boolean {
    return this.library.remove(id);
  }
  
  /**
   * Create fragment from molecule state
   */
  createFragment(
    state: MoleculeState,
    name: string,
    attachmentPoints: string[] = [],
    tags: string[] = []
  ): Fragment {
    return {
      id: crypto.randomUUID(),
      name,
      state: {
        atoms: new Map(state.atoms),
        bonds: new Map(state.bonds),
      },
      attachmentPoints,
      tags,
    };
  }
  
  /**
   * Insert fragment into molecule at attachment point
   */
  insertFragment(
    targetState: MoleculeState,
    fragment: Fragment,
    attachmentAtomId: string
  ): { atoms: string[]; bonds: string[] } {
    const insertedAtoms: string[] = [];
    const insertedBonds: string[] = [];
    
    // Copy fragment atoms (offset positions)
    const attachmentAtom = targetState.atoms.get(attachmentAtomId);
    if (!attachmentAtom) {
      return { atoms: [], bonds: [] };
    }
    
    const offsetX = attachmentAtom.x;
    const offsetY = attachmentAtom.y;
    
    fragment.state.atoms.forEach((fragAtom) => {
      const newAtomId = crypto.randomUUID();
      targetState.atoms.set(newAtomId, {
        id: newAtomId,
        element: fragAtom.element,
        x: fragAtom.x + offsetX + 50, // Offset to avoid overlap
        y: fragAtom.y + offsetY + 50,
        z: fragAtom.z || 0,
        charge: fragAtom.charge,
      });
      insertedAtoms.push(newAtomId);
    });
    
    // Copy fragment bonds
    const atomMapping = new Map<string, string>();
    let fragAtomIdx = 0;
    fragment.state.atoms.forEach((fragAtom) => {
      atomMapping.set(fragAtom.id, insertedAtoms[fragAtomIdx]);
      fragAtomIdx++;
    });
    
    fragment.state.bonds.forEach((fragBond) => {
      const newBondId = crypto.randomUUID();
      const atom1 = atomMapping.get(fragBond.atoms[0]);
      const atom2 = atomMapping.get(fragBond.atoms[1]);
      
      if (atom1 && atom2) {
        targetState.bonds.set(newBondId, {
          id: newBondId,
          atoms: [atom1, atom2],
          order: fragBond.order,
        });
        insertedBonds.push(newBondId);
      }
    });
    
    // Connect to attachment point if specified
    if (fragment.attachmentPoints.length > 0 && insertedAtoms.length > 0) {
      const fragAttachmentAtom = fragment.state.atoms.get(fragment.attachmentPoints[0]);
      if (fragAttachmentAtom) {
        const fragAtomIdx = Array.from(fragment.state.atoms.keys()).indexOf(fragAttachmentAtom.id);
        if (fragAtomIdx >= 0 && fragAtomIdx < insertedAtoms.length) {
          const newBondId = crypto.randomUUID();
          targetState.bonds.set(newBondId, {
            id: newBondId,
            atoms: [attachmentAtomId, insertedAtoms[fragAtomIdx]],
            order: 1,
          });
          insertedBonds.push(newBondId);
        }
      }
    }
    
    return { atoms: insertedAtoms, bonds: insertedBonds };
  }
}

export const fragmentEngine = new FragmentEngine();

