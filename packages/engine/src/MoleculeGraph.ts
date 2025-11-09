import { Atom, Bond, Element } from "./types";
import { getBondLength } from "./types";

export class MoleculeGraph {
  atoms: Map<string, Atom> = new Map();
  bonds: Map<string, Bond> = new Map();
  private nextAtomId = 1;
  private nextBondId = 1;

  addAtom(atom: Omit<Atom, "id">): string {
    const id = `atom_${this.nextAtomId++}`;
    this.atoms.set(id, { ...atom, id });
    return id;
  }

  removeAtom(id: string): boolean {
    if (!this.atoms.has(id)) return false;

    // Remove all bonds connected to this atom
    const bondsToRemove: string[] = [];
    this.bonds.forEach((bond, bondId) => {
      if (bond.a1 === id || bond.a2 === id) {
        bondsToRemove.push(bondId);
      }
    });
    bondsToRemove.forEach((bondId) => this.bonds.delete(bondId));

    this.atoms.delete(id);
    return true;
  }

  addBond(a1: string, a2: string, order: number = 1): string | null {
    if (a1 === a2) return null;
    if (!this.atoms.has(a1) || !this.atoms.has(a2)) return null;

    // Check if bond already exists
    const existingBond = Array.from(this.bonds.values()).find(
      (b) => (b.a1 === a1 && b.a2 === a2) || (b.a1 === a2 && b.a2 === a1)
    );
    if (existingBond) return existingBond.id;

    const id = `bond_${this.nextBondId++}`;
    this.bonds.set(id, { id, a1, a2, order });
    return id;
  }

  removeBond(id: string): boolean {
    return this.bonds.delete(id);
  }

  removeBondBetween(a1: string, a2: string): boolean {
    const bond = Array.from(this.bonds.values()).find(
      (b) => (b.a1 === a1 && b.a2 === a2) || (b.a1 === a2 && b.a2 === a1)
    );
    if (bond) {
      return this.bonds.delete(bond.id);
    }
    return false;
  }

  getNeighbors(atomId: string): string[] {
    const neighbors: string[] = [];
    this.bonds.forEach((bond) => {
      if (bond.a1 === atomId) neighbors.push(bond.a2);
      if (bond.a2 === atomId) neighbors.push(bond.a1);
    });
    return neighbors;
  }

  getBondOrder(a1: string, a2: string): number {
    const bond = Array.from(this.bonds.values()).find(
      (b) => (b.a1 === a1 && b.a2 === a2) || (b.a1 === a2 && b.a2 === a1)
    );
    return bond?.order || 0;
  }

  getBondsForAtom(atomId: string): Bond[] {
    return Array.from(this.bonds.values()).filter(
      (b) => b.a1 === atomId || b.a2 === atomId
    );
  }

  clone(): MoleculeGraph {
    const clone = new MoleculeGraph();
    this.atoms.forEach((atom) => {
      clone.atoms.set(atom.id, { ...atom });
    });
    this.bonds.forEach((bond) => {
      clone.bonds.set(bond.id, { ...bond });
    });
    clone.nextAtomId = this.nextAtomId;
    clone.nextBondId = this.nextBondId;
    return clone;
  }

  toJSON() {
    return {
      atoms: Array.from(this.atoms.values()),
      bonds: Array.from(this.bonds.values()),
    };
  }

  static fromJSON(data: { atoms: Atom[]; bonds: Bond[] }): MoleculeGraph {
    const graph = new MoleculeGraph();
    data.atoms.forEach((atom) => {
      graph.atoms.set(atom.id, atom);
    });
    data.bonds.forEach((bond) => {
      graph.bonds.set(bond.id, bond);
    });
    return graph;
  }

  getFormula(): string {
    const counts: Record<string, number> = {};
    this.atoms.forEach((atom) => {
      counts[atom.element] = (counts[atom.element] || 0) + 1;
    });

    // Standard order: C, H, then alphabetical
    const order = ["C", "H"];
    const other = Object.keys(counts)
      .filter((el) => !order.includes(el))
      .sort();

    return [...order, ...other]
      .filter((el) => counts[el])
      .map((el) => {
        const count = counts[el];
        return count === 1 ? el : `${el}${count}`;
      })
      .join("");
  }

  getMolecularWeight(): number {
    const atomicWeights: Record<string, number> = {
      H: 1.008,
      C: 12.011,
      N: 14.007,
      O: 15.999,
      F: 18.998,
      S: 32.065,
      P: 30.974,
      Cl: 35.453,
      Br: 79.904,
      I: 126.904,
    };

    let weight = 0;
    this.atoms.forEach((atom) => {
      weight += atomicWeights[atom.element] || 0;
    });
    return weight;
  }

  computeDegree(atomId: string): number {
    return this.getNeighbors(atomId).length;
  }
}

