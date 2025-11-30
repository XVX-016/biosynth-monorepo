/**
 * MoleculeStateEngine - Core molecule state management
 * 
 * This is the single source of truth for molecule structure.
 * Handles atoms, bonds, serialization, and basic validation.
 */

export interface Atom {
  id: string
  element: string
  x: number
  y: number
  z?: number
  charge: number
}

export interface Bond {
  id: string
  atoms: [string, string] // Atom IDs
  order: number
}

export class MoleculeStateEngine {
  atoms = new Map<string, Atom>()
  bonds = new Map<string, Bond>()

  /**
   * Add a new atom to the molecule
   */
  addAtom(element: string, x: number, y: number, z = 0, charge = 0): string {
    const id = crypto.randomUUID()
    this.atoms.set(id, { id, element, x, y, z, charge })
    return id
  }

  /**
   * Remove an atom and all its bonds
   */
  removeAtom(id: string): void {
    // Remove all bonds connected to this atom
    const bondsToRemove: string[] = []
    this.bonds.forEach((bond, bondId) => {
      if (bond.atoms[0] === id || bond.atoms[1] === id) {
        bondsToRemove.push(bondId)
      }
    })
    bondsToRemove.forEach(bondId => this.bonds.delete(bondId))
    
    // Remove the atom
    this.atoms.delete(id)
  }

  /**
   * Add a bond between two atoms
   */
  addBond(atomA: string, atomB: string, order = 1): string | null {
    // Validate atoms exist
    if (!this.atoms.has(atomA) || !this.atoms.has(atomB)) {
      return null
    }

    // Check if bond already exists
    for (const bond of this.bonds.values()) {
      if (
        (bond.atoms[0] === atomA && bond.atoms[1] === atomB) ||
        (bond.atoms[0] === atomB && bond.atoms[1] === atomA)
      ) {
        return null // Bond already exists
      }
    }

    const id = crypto.randomUUID()
    this.bonds.set(id, { id, atoms: [atomA, atomB], order })
    return id
  }

  /**
   * Remove a bond
   */
  removeBond(id: string): void {
    this.bonds.delete(id)
  }

  /**
   * Update atom position
   */
  setAtomPosition(id: string, x: number, y: number, z?: number): void {
    const atom = this.atoms.get(id)
    if (atom) {
      atom.x = x
      atom.y = y
      if (z !== undefined) atom.z = z
    }
  }

  /**
   * Get atom by ID
   */
  getAtom(id: string): Atom | undefined {
    return this.atoms.get(id)
  }

  /**
   * Get all atoms
   */
  getAllAtoms(): Atom[] {
    return Array.from(this.atoms.values())
  }

  /**
   * Get all bonds
   */
  getAllBonds(): Bond[] {
    return Array.from(this.bonds.values())
  }

  /**
   * Get bonds connected to an atom
   */
  getBondsForAtom(atomId: string): Bond[] {
    return Array.from(this.bonds.values()).filter(
      bond => bond.atoms[0] === atomId || bond.atoms[1] === atomId
    )
  }

  /**
   * Clear all atoms and bonds
   */
  clear(): void {
    this.atoms.clear()
    this.bonds.clear()
  }

  /**
   * Serialize to SMILES (requires backend)
   * TODO: Implement real SMILES generation via backend RDKit
   */
  async serializeSMILES(): Promise<string> {
    // For now, return placeholder
    // TODO: Call backend API to convert structure to SMILES
    return "C"
  }

  /**
   * Serialize to SDF format
   * TODO: Implement real SDF serialization
   */
  toSDF(): string {
    // Basic SDF structure
    let sdf = ""
    
    // Header
    sdf += "\n\n"
    
    // Counts line: atoms bonds
    sdf += `${this.atoms.size.toString().padStart(3, '0')}${this.bonds.size.toString().padStart(3, '0')}  0  0  0  0  0  0  0  0  1 V2000\n`
    
    // Atom block
    const atomArray = Array.from(this.atoms.values())
    atomArray.forEach(atom => {
      const x = atom.x.toFixed(4).padStart(10)
      const y = atom.y.toFixed(4).padStart(10)
      const z = (atom.z || 0).toFixed(4).padStart(10)
      const element = atom.element.padEnd(3)
      sdf += `${x}${y}${z} ${element} 0  0  0  0  0  0  0  0  0  0  0  0\n`
    })
    
    // Bond block
    this.bonds.forEach(bond => {
      const idxA = atomArray.findIndex(a => a.id === bond.atoms[0])
      const idxB = atomArray.findIndex(a => a.id === bond.atoms[1])
      if (idxA >= 0 && idxB >= 0) {
        const a = (idxA + 1).toString().padStart(3)
        const b = (idxB + 1).toString().padStart(3)
        const order = bond.order.toString().padStart(3)
        sdf += `${a}${b}${order}  0  0  0  0\n`
      }
    })
    
    sdf += "M  END\n$$$$\n"
    return sdf
  }

  /**
   * Load from SMILES (requires backend)
   * TODO: Implement real SMILES parsing via backend
   */
  async loadFromSMILES(smiles: string): Promise<void> {
    // TODO: Call backend API to parse SMILES and populate atoms/bonds
    this.clear()
  }

  /**
   * Export to JSON for persistence
   */
  toJSON(): { atoms: Atom[], bonds: Bond[] } {
    return {
      atoms: this.getAllAtoms(),
      bonds: this.getAllBonds()
    }
  }

  /**
   * Load from JSON
   */
  fromJSON(data: { atoms: Atom[], bonds: Bond[] }): void {
    this.clear()
    data.atoms.forEach(atom => {
      this.atoms.set(atom.id, atom)
    })
    data.bonds.forEach(bond => {
      this.bonds.set(bond.id, bond)
    })
  }
}

// Singleton instance
export const moleculeEngine = new MoleculeStateEngine()

