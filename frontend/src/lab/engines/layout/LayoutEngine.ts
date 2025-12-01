/**
 * LayoutEngine - Generate 2D coordinates for molecules
 */

import type { MoleculeState, Atom } from '../MoleculeStateEngine';

export interface LayoutOptions {
  spacing?: number;
  center?: { x: number; y: number };
  algorithm?: 'force-directed' | 'ring-based' | 'tree';
}

/**
 * Simple force-directed layout algorithm
 */
function forceDirectedLayout(
  state: MoleculeState,
  options: LayoutOptions
): Map<string, { x: number; y: number }> {
  const spacing = options.spacing ?? 50;
  const positions = new Map<string, { x: number; y: number }>();
  const atoms = Array.from(state.atoms.values());
  const bonds = Array.from(state.bonds.values());
  
  // Initialize positions randomly or from existing
  atoms.forEach((atom, idx) => {
    if (atom.x === 0 && atom.y === 0) {
      // Place in a circle if no position
      const angle = (idx / atoms.length) * Math.PI * 2;
      positions.set(atom.id, {
        x: Math.cos(angle) * spacing * 2,
        y: Math.sin(angle) * spacing * 2,
      });
    } else {
      positions.set(atom.id, { x: atom.x, y: atom.y });
    }
  });
  
  // Simple force-directed iterations
  const iterations = 100;
  const k = spacing; // Spring constant
  const repulsion = spacing * 2;
  
  for (let iter = 0; iter < iterations; iter++) {
    const forces = new Map<string, { fx: number; fy: number }>();
    
    // Initialize forces
    atoms.forEach((atom) => {
      forces.set(atom.id, { fx: 0, fy: 0 });
    });
    
    // Repulsion between all atoms
    for (let i = 0; i < atoms.length; i++) {
      for (let j = i + 1; j < atoms.length; j++) {
        const a1 = atoms[i];
        const a2 = atoms[j];
        const pos1 = positions.get(a1.id)!;
        const pos2 = positions.get(a2.id)!;
        
        const dx = pos2.x - pos1.x;
        const dy = pos2.y - pos1.y;
        const dist = Math.sqrt(dx * dx + dy * dy) || 1;
        
        const force = repulsion / (dist * dist);
        const fx = (dx / dist) * force;
        const fy = (dy / dist) * force;
        
        const f1 = forces.get(a1.id)!;
        const f2 = forces.get(a2.id)!;
        f1.fx -= fx;
        f1.fy -= fy;
        f2.fx += fx;
        f2.fy += fy;
      }
    }
    
    // Attraction along bonds
    bonds.forEach((bond) => {
      const pos1 = positions.get(bond.atoms[0])!;
      const pos2 = positions.get(bond.atoms[1])!;
      
      const dx = pos2.x - pos1.x;
      const dy = pos2.y - pos1.y;
      const dist = Math.sqrt(dx * dx + dy * dy) || 1;
      
      const force = (dist - spacing) * k * 0.01;
      const fx = (dx / dist) * force;
      const fy = (dy / dist) * force;
      
      const f1 = forces.get(bond.atoms[0])!;
      const f2 = forces.get(bond.atoms[1])!;
      f1.fx += fx;
      f1.fy += fy;
      f2.fx -= fx;
      f2.fy -= fy;
    });
    
    // Apply forces
    atoms.forEach((atom) => {
      const pos = positions.get(atom.id)!;
      const force = forces.get(atom.id)!;
      
      pos.x += force.fx * 0.1;
      pos.y += force.fy * 0.1;
    });
  }
  
  // Center the molecule
  if (options.center) {
    const centerX = Array.from(positions.values()).reduce((sum, p) => sum + p.x, 0) / positions.size;
    const centerY = Array.from(positions.values()).reduce((sum, p) => sum + p.y, 0) / positions.size;
    
    positions.forEach((pos) => {
      pos.x += options.center!.x - centerX;
      pos.y += options.center!.y - centerY;
    });
  }
  
  return positions;
}

/**
 * Ring-based layout for cyclic molecules
 */
function ringBasedLayout(
  state: MoleculeState,
  options: LayoutOptions
): Map<string, { x: number; y: number }> {
  const spacing = options.spacing ?? 50;
  const positions = new Map<string, { x: number; y: number }>();
  
  // Detect rings
  const rings = detectRings(state);
  
  if (rings.length > 0) {
    // Place first ring in center
    const firstRing = rings[0];
    const angleStep = (Math.PI * 2) / firstRing.length;
    
    firstRing.forEach((atomId, idx) => {
      const angle = idx * angleStep;
      positions.set(atomId, {
        x: Math.cos(angle) * spacing,
        y: Math.sin(angle) * spacing,
      });
    });
    
    // Place other atoms using force-directed
    return forceDirectedLayout(state, options);
  }
  
  return forceDirectedLayout(state, options);
}

/**
 * Detect rings in molecule
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
  
  // Remove duplicates
  return rings.filter((ring, index, self) => {
    const normalized = ring.sort().join('-');
    return index === self.findIndex((r) => r.sort().join('-') === normalized);
  });
}

/**
 * Main LayoutEngine class
 */
export class LayoutEngine {
  /**
   * Generate 2D coordinates for molecule
   */
  generate2DCoordinates(
    state: MoleculeState,
    options: LayoutOptions = {}
  ): Map<string, { x: number; y: number }> {
    const algorithm = options.algorithm || 'force-directed';
    
    switch (algorithm) {
      case 'ring-based':
        return ringBasedLayout(state, options);
      case 'tree':
        // TODO: Implement tree layout
        return forceDirectedLayout(state, options);
      case 'force-directed':
      default:
        return forceDirectedLayout(state, options);
    }
  }
  
  /**
   * Apply layout to molecule state (mutates state)
   */
  applyLayout(state: MoleculeState, options?: LayoutOptions): void {
    const positions = this.generate2DCoordinates(state, options);
    
    positions.forEach((pos, atomId) => {
      const atom = state.atoms.get(atomId);
      if (atom) {
        atom.x = pos.x;
        atom.y = pos.y;
      }
    });
  }
}

export const layoutEngine = new LayoutEngine();

