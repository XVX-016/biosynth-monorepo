/**
 * Geometry3D - Generate 3D coordinates for molecules
 */

import type { MoleculeState, Atom } from '../MoleculeStateEngine';

export interface Geometry3DOptions {
  method?: 'force-field' | 'distance-geometry' | 'template';
  optimize?: boolean;
}

/**
 * Simple distance geometry for 3D coordinates
 */
function distanceGeometry3D(
  state: MoleculeState,
  options: Geometry3DOptions
): Map<string, { x: number; y: number; z: number }> {
  const positions = new Map<string, { x: number; y: number; z: number }>();
  const atoms = Array.from(state.atoms.values());
  const bonds = Array.from(state.bonds.values());
  
  // Bond length table (in Angstroms)
  const BOND_LENGTHS: Record<string, number> = {
    'C-C': 1.54,
    'C=C': 1.34,
    'C≡C': 1.20,
    'C-H': 1.09,
    'C-O': 1.43,
    'C=O': 1.23,
    'C-N': 1.47,
    'C=N': 1.29,
    'O-H': 0.96,
    'N-H': 1.01,
    'N-N': 1.45,
    'N=N': 1.25,
  };
  
  function getBondLength(element1: string, element2: string, order: number): number {
    const key = order === 1 ? `${element1}-${element2}` :
                order === 2 ? `${element1}=${element2}` :
                `${element1}≡${element2}`;
    return BOND_LENGTHS[key] || BOND_LENGTHS[`${element2}-${element1}`] || 1.5;
  }
  
  // Initialize positions from 2D or randomly
  atoms.forEach((atom, idx) => {
    if (atom.z === undefined || atom.z === 0) {
      // Generate 3D position from 2D
      positions.set(atom.id, {
        x: atom.x / 20, // Scale down
        y: atom.y / 20,
        z: (Math.random() - 0.5) * 0.5, // Small random z
      });
    } else {
      positions.set(atom.id, {
        x: atom.x / 20,
        y: atom.y / 20,
        z: atom.z / 20,
      });
    }
  });
  
  // Distance geometry: adjust positions to match bond lengths
  const iterations = 50;
  
  for (let iter = 0; iter < iterations; iter++) {
    bonds.forEach((bond) => {
      const atom1 = state.atoms.get(bond.atoms[0]);
      const atom2 = state.atoms.get(bond.atoms[1]);
      
      if (!atom1 || !atom2) return;
      
      const pos1 = positions.get(atom1.id)!;
      const pos2 = positions.get(atom2.id)!;
      
      const dx = pos2.x - pos1.x;
      const dy = pos2.y - pos1.y;
      const dz = pos2.z - pos1.z;
      const currentDist = Math.sqrt(dx * dx + dy * dy + dz * dz);
      
      const targetDist = getBondLength(atom1.element, atom2.element, bond.order);
      const diff = targetDist - currentDist;
      
      if (Math.abs(diff) > 0.01) {
        const scale = diff / (currentDist || 1);
        const moveX = dx * scale * 0.5;
        const moveY = dy * scale * 0.5;
        const moveZ = dz * scale * 0.5;
        
        pos1.x -= moveX;
        pos1.y -= moveY;
        pos1.z -= moveZ;
        pos2.x += moveX;
        pos2.y += moveY;
        pos2.z += moveZ;
      }
    });
    
    // Repulsion between non-bonded atoms
    for (let i = 0; i < atoms.length; i++) {
      for (let j = i + 1; j < atoms.length; j++) {
        const a1 = atoms[i];
        const a2 = atoms[j];
        
        // Check if bonded
        const isBonded = bonds.some(
          (bond) =>
            (bond.atoms[0] === a1.id && bond.atoms[1] === a2.id) ||
            (bond.atoms[0] === a2.id && bond.atoms[1] === a1.id)
        );
        
        if (!isBonded) {
          const pos1 = positions.get(a1.id)!;
          const pos2 = positions.get(a2.id)!;
          
          const dx = pos2.x - pos1.x;
          const dy = pos2.y - pos1.y;
          const dz = pos2.z - pos1.z;
          const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
          
          const minDist = 2.0; // Minimum distance between non-bonded atoms
          if (dist < minDist && dist > 0) {
            const force = (minDist - dist) * 0.1;
            const scale = force / (dist || 1);
            
            pos1.x -= dx * scale * 0.5;
            pos1.y -= dy * scale * 0.5;
            pos1.z -= dz * scale * 0.5;
            pos2.x += dx * scale * 0.5;
            pos2.y += dy * scale * 0.5;
            pos2.z += dz * scale * 0.5;
          }
        }
      }
    }
  }
  
  return positions;
}

/**
 * Main Geometry3D class
 */
export class Geometry3D {
  /**
   * Generate 3D coordinates
   */
  generate3D(
    state: MoleculeState,
    options: Geometry3DOptions = {}
  ): Map<string, { x: number; y: number; z: number }> {
    const method = options.method || 'distance-geometry';
    
    switch (method) {
      case 'distance-geometry':
        return distanceGeometry3D(state, options);
      case 'force-field':
        // TODO: Implement force-field optimization
        return distanceGeometry3D(state, options);
      case 'template':
        // TODO: Implement template-based geometry
        return distanceGeometry3D(state, options);
      default:
        return distanceGeometry3D(state, options);
    }
  }
  
  /**
   * Apply 3D geometry to molecule state
   */
  apply3D(state: MoleculeState, options?: Geometry3DOptions): void {
    const positions = this.generate3D(state, options);
    
    positions.forEach((pos, atomId) => {
      const atom = state.atoms.get(atomId);
      if (atom) {
        atom.x = pos.x * 20; // Scale up for display
        atom.y = pos.y * 20;
        atom.z = pos.z * 20;
      }
    });
  }
}

export const geometry3D = new Geometry3D();

