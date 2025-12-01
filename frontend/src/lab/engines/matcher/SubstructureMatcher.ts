/**
 * SubstructureMatcher - Match SMARTS patterns in molecules
 */

import type { MoleculeState, Atom, Bond } from '../MoleculeStateEngine';

export interface MatchResult {
  matches: Array<{
    atoms: string[];
    bonds: string[];
    score: number;
  }>;
  pattern: string;
}

/**
 * Simple SMARTS pattern parser (basic implementation)
 * Full SMARTS requires backend RDKit
 */
class SMARTSParser {
  private pattern: string;
  
  constructor(pattern: string) {
    this.pattern = pattern;
  }
  
  /**
   * Parse basic SMARTS patterns
   * Supports: [C], [N], [O], [*] (any), bonds: -, =, #
   */
  parse(): any {
    // TODO: Full SMARTS parsing requires RDKit
    // This is a simplified version for common patterns
    
    const atoms: string[] = [];
    const bonds: number[] = [];
    
    let i = 0;
    while (i < this.pattern.length) {
      if (this.pattern[i] === '[') {
        // Atom pattern
        const end = this.pattern.indexOf(']', i);
        if (end > i) {
          const atomPattern = this.pattern.substring(i + 1, end);
          atoms.push(atomPattern);
          i = end + 1;
        } else {
          i++;
        }
      } else if (this.pattern[i] === '-') {
        bonds.push(1);
        i++;
      } else if (this.pattern[i] === '=') {
        bonds.push(2);
        i++;
      } else if (this.pattern[i] === '#') {
        bonds.push(3);
        i++;
      } else if (this.pattern[i] === '*') {
        atoms.push('*'); // Any atom
        i++;
      } else {
        // Single character element
        atoms.push(this.pattern[i]);
        i++;
      }
    }
    
    return { atoms, bonds };
  }
}

/**
 * Match atom pattern
 */
function matchAtom(atom: Atom, pattern: string): boolean {
  if (pattern === '*' || pattern === '[*]') return true;
  if (pattern === atom.element) return true;
  
  // Element groups
  if (pattern === '[C]' && atom.element === 'C') return true;
  if (pattern === '[N]' && atom.element === 'N') return true;
  if (pattern === '[O]' && atom.element === 'O') return true;
  if (pattern === '[H]' && atom.element === 'H') return true;
  
  return false;
}

/**
 * Match bond pattern
 */
function matchBond(bond: Bond, pattern: number): boolean {
  return bond.order === pattern;
}

/**
 * Find subgraph isomorphism (simplified)
 */
function findSubgraph(
  state: MoleculeState,
  patternAtoms: string[],
  patternBonds: number[]
): MatchResult['matches'] {
  const matches: MatchResult['matches'] = [];
  const atoms = Array.from(state.atoms.values());
  const bonds = Array.from(state.bonds.values());
  
  // Simple matching: find sequences of atoms matching pattern
  // Full implementation would use VF2 or Ullmann algorithm
  
  if (patternAtoms.length === 1) {
    // Single atom pattern
    atoms.forEach((atom) => {
      if (matchAtom(atom, patternAtoms[0])) {
        matches.push({
          atoms: [atom.id],
          bonds: [],
          score: 1.0,
        });
      }
    });
  } else if (patternAtoms.length === 2) {
    // Two-atom pattern with bond
    const bondOrder = patternBonds[0] || 1;
    
    bonds.forEach((bond) => {
      const atom1 = state.atoms.get(bond.atoms[0]);
      const atom2 = state.atoms.get(bond.atoms[1]);
      
      if (!atom1 || !atom2) return;
      
      if (
        matchAtom(atom1, patternAtoms[0]) &&
        matchAtom(atom2, patternAtoms[1]) &&
        matchBond(bond, bondOrder)
      ) {
        matches.push({
          atoms: [atom1.id, atom2.id],
          bonds: [bond.id],
          score: 1.0,
        });
      }
      
      // Try reverse
      if (
        matchAtom(atom1, patternAtoms[1]) &&
        matchAtom(atom2, patternAtoms[0]) &&
        matchBond(bond, bondOrder)
      ) {
        matches.push({
          atoms: [atom2.id, atom1.id],
          bonds: [bond.id],
          score: 1.0,
        });
      }
    });
  } else {
    // Multi-atom pattern - use DFS
    const visited = new Set<string>();
    
    const dfs = (
      currentAtom: Atom,
      patternIdx: number,
      matchedAtoms: string[],
      matchedBonds: string[]
    ) => {
      if (patternIdx >= patternAtoms.length) {
        // Found a match
        matches.push({
          atoms: [...matchedAtoms],
          bonds: [...matchedBonds],
          score: 1.0,
        });
        return;
      }
      
      if (!matchAtom(currentAtom, patternAtoms[patternIdx])) {
        return;
      }
      
      matchedAtoms.push(currentAtom.id);
      visited.add(currentAtom.id);
      
      if (patternIdx < patternAtoms.length - 1) {
        // Need to find next atom
        const nextBondOrder = patternBonds[patternIdx] || 1;
        const connectedBonds = bonds.filter(
          (bond) =>
            (bond.atoms[0] === currentAtom.id || bond.atoms[1] === currentAtom.id) &&
            bond.order === nextBondOrder &&
            !matchedBonds.includes(bond.id)
        );
        
        connectedBonds.forEach((bond) => {
          const nextAtomId = bond.atoms[0] === currentAtom.id ? bond.atoms[1] : bond.atoms[0];
          const nextAtom = state.atoms.get(nextAtomId);
          
          if (nextAtom && !visited.has(nextAtomId)) {
            matchedBonds.push(bond.id);
            dfs(nextAtom, patternIdx + 1, matchedAtoms, matchedBonds);
            matchedBonds.pop();
          }
        });
      } else {
        // Last atom in pattern
        matches.push({
          atoms: [...matchedAtoms],
          bonds: [...matchedBonds],
          score: 1.0,
        });
      }
      
      matchedAtoms.pop();
      visited.delete(currentAtom.id);
    };
    
    // Try starting from each atom
    atoms.forEach((atom) => {
      dfs(atom, 0, [], []);
    });
  }
  
  // Remove duplicates
  const uniqueMatches = matches.filter((match, index, self) => {
    const key = match.atoms.sort().join('-');
    return index === self.findIndex((m) => m.atoms.sort().join('-') === key);
  });
  
  return uniqueMatches;
}

/**
 * Main SubstructureMatcher class
 */
export class SubstructureMatcher {
  /**
   * Match SMARTS pattern in molecule
   */
  matchSMARTS(state: MoleculeState, pattern: string): MatchResult {
    const parser = new SMARTSParser(pattern);
    const parsed = parser.parse();
    
    const matches = findSubgraph(state, parsed.atoms, parsed.bonds);
    
    return {
      matches,
      pattern,
    };
  }
  
  /**
   * Check if pattern exists in molecule
   */
  hasPattern(state: MoleculeState, pattern: string): boolean {
    const result = this.matchSMARTS(state, pattern);
    return result.matches.length > 0;
  }
  
  /**
   * Count occurrences of pattern
   */
  countPattern(state: MoleculeState, pattern: string): number {
    const result = this.matchSMARTS(state, pattern);
    return result.matches.length;
  }
}

export const substructureMatcher = new SubstructureMatcher();

