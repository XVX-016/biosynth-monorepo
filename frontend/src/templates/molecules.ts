// src/templates/molecules.ts
import { MoleculeGraph } from "@biosynth/engine";
import { placeAtomRelative } from "@biosynth/engine";

/**
 * Create a benzene ring (C6H6)
 */
export function createBenzene(): MoleculeGraph {
  const g = new MoleculeGraph();
  const atomIds: string[] = [];

  // Create 6 carbon atoms in a ring
  const center = [0, 0, 0];
  const radius = 1.4; // Approximate benzene ring radius
  
  for (let i = 0; i < 6; i++) {
    const angle = (i * Math.PI * 2) / 6;
    const x = center[0] + radius * Math.cos(angle);
    const y = center[1] + radius * Math.sin(angle);
    const z = center[2];
    
    const id = g.addAtom({ element: "C", position: [x, y, z] });
    atomIds.push(id);
  }

  // Create alternating single/double bonds (Kekulé structure)
  for (let i = 0; i < 6; i++) {
    const a = atomIds[i];
    const b = atomIds[(i + 1) % 6];
    const order = i % 2 === 0 ? 2 : 1; // Alternating double/single
    g.addBond(a, b, order);
  }

  // Add hydrogen atoms
  const hAtomIds: string[] = [];
  for (let i = 0; i < 6; i++) {
    const cAtom = g.atoms.get(atomIds[i])!;
    const hPosition = placeAtomRelative(cAtom, "H");
    const hId = g.addAtom({ element: "H", position: hPosition });
    hAtomIds.push(hId);
    g.addBond(atomIds[i], hId, 1);
  }

  return g;
}

/**
 * Create methane (CH4)
 */
export function createMethane(): MoleculeGraph {
  const g = new MoleculeGraph();
  
  const cId = g.addAtom({ element: "C", position: [0, 0, 0] });
  const cAtom = g.atoms.get(cId)!;

  // Tetrahedral arrangement
  const h1Id = g.addAtom({ element: "H", position: placeAtomRelative(cAtom, "H", [1, 1, 1]) });
  const h2Id = g.addAtom({ element: "H", position: placeAtomRelative(cAtom, "H", [-1, -1, 1]) });
  const h3Id = g.addAtom({ element: "H", position: placeAtomRelative(cAtom, "H", [1, -1, -1]) });
  const h4Id = g.addAtom({ element: "H", position: placeAtomRelative(cAtom, "H", [-1, 1, -1]) });

  g.addBond(cId, h1Id, 1);
  g.addBond(cId, h2Id, 1);
  g.addBond(cId, h3Id, 1);
  g.addBond(cId, h4Id, 1);

  return g;
}

/**
 * Create water (H2O)
 */
export function createWater(): MoleculeGraph {
  const g = new MoleculeGraph();
  
  const oId = g.addAtom({ element: "O", position: [0, 0, 0] });
  const oAtom = g.atoms.get(oId)!;

  const h1Id = g.addAtom({ element: "H", position: placeAtomRelative(oAtom, "H", [0.96, 0.3, 0]) });
  const h2Id = g.addAtom({ element: "H", position: placeAtomRelative(oAtom, "H", [-0.96, 0.3, 0]) });

  g.addBond(oId, h1Id, 1);
  g.addBond(oId, h2Id, 1);

  return g;
}

/**
 * Create ethanol (C2H5OH)
 */
export function createEthanol(): MoleculeGraph {
  const g = new MoleculeGraph();
  
  // First carbon
  const c1Id = g.addAtom({ element: "C", position: [-0.7, 0, 0] });
  const c1Atom = g.atoms.get(c1Id)!;

  // Second carbon
  const c2Id = g.addAtom({ element: "C", position: [0.7, 0, 0] });
  const c2Atom = g.atoms.get(c2Id)!;

  // Oxygen
  const oId = g.addAtom({ element: "O", position: placeAtomRelative(c2Atom, "O", [0, 0.8, 0]) });
  const oAtom = g.atoms.get(oId)!;

  // Hydrogens on C1
  const h1Id = g.addAtom({ element: "H", position: placeAtomRelative(c1Atom, "H", [-1.2, -0.9, 0]) });
  const h2Id = g.addAtom({ element: "H", position: placeAtomRelative(c1Atom, "H", [-1.2, 0.9, 0]) });
  const h3Id = g.addAtom({ element: "H", position: placeAtomRelative(c1Atom, "H", [0, -0.9, 0]) });

  // Hydrogens on C2
  const h4Id = g.addAtom({ element: "H", position: placeAtomRelative(c2Atom, "H", [0, -0.9, 0]) });
  const h5Id = g.addAtom({ element: "H", position: placeAtomRelative(c2Atom, "H", [0, 0.9, 0]) });

  // Hydrogen on O
  const h6Id = g.addAtom({ element: "H", position: placeAtomRelative(oAtom, "H", [0, 1.4, 0]) });

  // Create bonds
  g.addBond(c1Id, c2Id, 1);
  g.addBond(c2Id, oId, 1);
  g.addBond(c1Id, h1Id, 1);
  g.addBond(c1Id, h2Id, 1);
  g.addBond(c1Id, h3Id, 1);
  g.addBond(c2Id, h4Id, 1);
  g.addBond(c2Id, h5Id, 1);
  g.addBond(oId, h6Id, 1);

  return g;
}

/**
 * Create acetic acid (CH3COOH)
 */
export function createAceticAcid(): MoleculeGraph {
  const g = new MoleculeGraph();
  
  // Methyl carbon
  const c1Id = g.addAtom({ element: "C", position: [-1.5, 0, 0] });
  const c1Atom = g.atoms.get(c1Id)!;

  // Carbonyl carbon
  const c2Id = g.addAtom({ element: "C", position: [0, 0, 0] });
  const c2Atom = g.atoms.get(c2Id)!;

  // Carbonyl oxygen
  const o1Id = g.addAtom({ element: "O", position: placeAtomRelative(c2Atom, "O", [0.7, 0.7, 0]) });

  // Hydroxyl oxygen
  const o2Id = g.addAtom({ element: "O", position: placeAtomRelative(c2Atom, "O", [-0.7, 0.7, 0]) });
  const o2Atom = g.atoms.get(o2Id)!;

  // Hydrogens on methyl
  const h1Id = g.addAtom({ element: "H", position: placeAtomRelative(c1Atom, "H", [-2.0, -0.9, 0]) });
  const h2Id = g.addAtom({ element: "H", position: placeAtomRelative(c1Atom, "H", [-2.0, 0.9, 0]) });
  const h3Id = g.addAtom({ element: "H", position: placeAtomRelative(c1Atom, "H", [-1.0, 0, 0]) });

  // Hydrogen on hydroxyl
  const h4Id = g.addAtom({ element: "H", position: placeAtomRelative(o2Atom, "H", [-1.2, 1.4, 0]) });

  // Create bonds
  g.addBond(c1Id, c2Id, 1);
  g.addBond(c2Id, o1Id, 2); // Double bond
  g.addBond(c2Id, o2Id, 1);
  g.addBond(c1Id, h1Id, 1);
  g.addBond(c1Id, h2Id, 1);
  g.addBond(c1Id, h3Id, 1);
  g.addBond(o2Id, h4Id, 1);

  return g;
}

