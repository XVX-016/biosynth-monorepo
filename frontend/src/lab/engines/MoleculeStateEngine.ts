/**
 * MoleculeStateEngine - Core molecule state management
 * 
 * This is the single source of truth for molecule data.
 * All editing operations go through this engine.
 */

export interface Atom {
  id: string;
  element: string;
  x: number;
  y: number;
  z?: number;
  charge: number;
}

export interface Bond {
  id: string;
  atoms: [string, string];
  order: number;
}

export interface MoleculeState {
  atoms: Map<string, Atom>;
  bonds: Map<string, Bond>;
}

export class MoleculeStateEngine {
  atoms = new Map<string, Atom>();
  bonds = new Map<string, Bond>();

  /**
   * Add an atom to the molecule
   */
  addAtom(element: string, x: number, y: number, z: number = 0, charge: number = 0): string {
    const id = crypto.randomUUID();
    this.atoms.set(id, { id, element, x, y, z, charge });
    return id;
  }

  /**
   * Remove an atom and all its bonds
   */
  removeAtom(atomId: string): void {
    // Remove all bonds connected to this atom
    const bondsToRemove: string[] = [];
    this.bonds.forEach((bond, bondId) => {
      if (bond.atoms[0] === atomId || bond.atoms[1] === atomId) {
        bondsToRemove.push(bondId);
      }
    });
    bondsToRemove.forEach(bondId => this.bonds.delete(bondId));
    
    // Remove the atom
    this.atoms.delete(atomId);
  }

  /**
   * Add a bond between two atoms
   */
  addBond(atom1Id: string, atom2Id: string, order: number = 1): string {
    if (atom1Id === atom2Id) {
      throw new Error('Cannot bond atom to itself');
    }
    
    // Check if bond already exists
    for (const bond of this.bonds.values()) {
      if (
        (bond.atoms[0] === atom1Id && bond.atoms[1] === atom2Id) ||
        (bond.atoms[0] === atom2Id && bond.atoms[1] === atom1Id)
      ) {
        throw new Error('Bond already exists');
      }
    }

    const id = crypto.randomUUID();
    this.bonds.set(id, { id, atoms: [atom1Id, atom2Id], order });
    return id;
  }

  /**
   * Remove a bond
   */
  removeBond(bondId: string): void {
    this.bonds.delete(bondId);
  }

  /**
   * Update atom position
   */
  setAtomPosition(atomId: string, x: number, y: number, z?: number): void {
    const atom = this.atoms.get(atomId);
    if (atom) {
      atom.x = x;
      atom.y = y;
      if (z !== undefined) atom.z = z;
    }
  }

  /**
   * Update atom element
   */
  setAtomElement(atomId: string, element: string): void {
    const atom = this.atoms.get(atomId);
    if (atom) {
      atom.element = element;
    }
  }

  /**
   * Update bond order
   */
  setBondOrder(bondId: string, order: number): void {
    const bond = this.bonds.get(bondId);
    if (bond) {
      bond.order = order;
    }
  }

  /**
   * Get all atoms as array
   */
  getAtoms(): Atom[] {
    return Array.from(this.atoms.values());
  }

  /**
   * Get all bonds as array
   */
  getBonds(): Bond[] {
    return Array.from(this.bonds.values());
  }

  /**
   * Get atom by ID
   */
  getAtom(atomId: string): Atom | undefined {
    return this.atoms.get(atomId);
  }

  /**
   * Get bond by ID
   */
  getBond(bondId: string): Bond | undefined {
    return this.bonds.get(bondId);
  }

  /**
   * Get bonds connected to an atom
   */
  getAtomBonds(atomId: string): Bond[] {
    return this.getBonds().filter(
      bond => bond.atoms[0] === atomId || bond.atoms[1] === atomId
    );
  }

  /**
   * Clear all atoms and bonds
   */
  clear(): void {
    this.atoms.clear();
    this.bonds.clear();
  }

  /**
   * Get current state snapshot
   */
  getState(): MoleculeState {
    return {
      atoms: new Map(this.atoms),
      bonds: new Map(this.bonds),
    };
  }

  /**
   * Load state from snapshot
   */
  loadState(state: MoleculeState): void {
    this.atoms = new Map(state.atoms);
    this.bonds = new Map(state.bonds);
  }

  /**
   * Serialize to SMILES (placeholder - will call backend)
   */
  serializeSMILES(): string {
    // TODO: Real implementation via backend RDKit
    return "C";
  }

  /**
   * Serialize to SDF format (placeholder)
   */
  toSDF(): string {
    // TODO: Implement real SDF serialization
    let sdf = "";
    
    // Header
    sdf += "\n\n\n";
    
    // Counts line
    sdf += `  ${this.atoms.size.toString().padStart(3, '0')}  ${this.bonds.size.toString().padStart(3, '0')}  0  0  0  0  0  0  0  0999 V2000\n`;
    
    // Atom block
    const atomArray = Array.from(this.atoms.values());
    atomArray.forEach(atom => {
      const x = atom.x.toFixed(4).padStart(10);
      const y = atom.y.toFixed(4).padStart(10);
      const z = (atom.z || 0).toFixed(4).padStart(10);
      const element = atom.element.padEnd(3);
      sdf += `${x}${y}${z} ${element} 0  0  0  0  0  0  0  0  0  0\n`;
    });
    
    // Bond block
    this.bonds.forEach(bond => {
      const idx1 = atomArray.findIndex(a => a.id === bond.atoms[0]) + 1;
      const idx2 = atomArray.findIndex(a => a.id === bond.atoms[1]) + 1;
      sdf += `${idx1.toString().padStart(3)}${idx2.toString().padStart(3)}${bond.order}  0  0  0  0\n`;
    });
    
    sdf += "M  END\n$$$$\n";
    return sdf;
  }
}

export const moleculeEngine = new MoleculeStateEngine();

