// src/templates/groups.ts
import { MoleculeGraph } from "@biosynth/engine";
import { atomTemplates } from "./atoms";
import { placeAtomRelative } from "@biosynth/engine";

export interface TemplateInstance {
  atomIds: string[]; 
  entryAtomId: string; // The atom that connects to the parent molecule
}

/**
 * Apply a hydroxyl group (-OH) to an atom
 */
export function applyHydroxyl(graph: MoleculeGraph, attachToId: string): TemplateInstance | null {
  const attachAtom = graph.atoms.get(attachToId);
  if (!attachAtom) return null;

  // Place O relative to attachment point
  const oPosition = placeAtomRelative(attachAtom, "O");
  const oId = graph.addAtom({ element: "O", position: oPosition });

  // Place H relative to O
  const oAtom = graph.atoms.get(oId)!;
  const hPosition = placeAtomRelative(oAtom, "H");
  const hId = graph.addAtom({ element: "H", position: hPosition });

  // Create bonds (use addBond directly since we control placement)
  graph.addBond(attachToId, oId, 1);
  graph.addBond(oId, hId, 1);

  return { atomIds: [oId, hId], entryAtomId: oId };
}

/**
 * Apply a methyl group (-CH3) to an atom
 */
export function applyMethyl(graph: MoleculeGraph, attachToId: string): TemplateInstance | null {
  const attachAtom = graph.atoms.get(attachToId);
  if (!attachAtom) return null;

  // Place C relative to attachment point
  const cPosition = placeAtomRelative(attachAtom, "C");
  const cId = graph.addAtom({ element: "C", position: cPosition });
  const cAtom = graph.atoms.get(cId)!;

  // Place three H atoms around C (tetrahedral-like)
  const h1Position = placeAtomRelative(cAtom, "H", [1, 1, 0]);
  const h1Id = graph.addAtom({ element: "H", position: h1Position });

  const h2Position = placeAtomRelative(cAtom, "H", [-1, 1, 0]);
  const h2Id = graph.addAtom({ element: "H", position: h2Position });

  const h3Position = placeAtomRelative(cAtom, "H", [0, -1, 1]);
  const h3Id = graph.addAtom({ element: "H", position: h3Position });

  // Create bonds (use addBond directly since we control placement)
  graph.addBond(attachToId, cId, 1);
  graph.addBond(cId, h1Id, 1);
  graph.addBond(cId, h2Id, 1);
  graph.addBond(cId, h3Id, 1);

  return { atomIds: [cId, h1Id, h2Id, h3Id], entryAtomId: cId };
}

/**
 * Apply a carbonyl group (C=O) to an atom
 */
export function applyCarbonyl(graph: MoleculeGraph, attachToId: string): TemplateInstance | null {
  const attachAtom = graph.atoms.get(attachToId);
  if (!attachAtom) return null;

  // Place C relative to attachment point
  const cPosition = placeAtomRelative(attachAtom, "C");
  const cId = graph.addAtom({ element: "C", position: cPosition });
  const cAtom = graph.atoms.get(cId)!;

  // Place O relative to C
  const oPosition = placeAtomRelative(cAtom, "O");
  const oId = graph.addAtom({ element: "O", position: oPosition });

  // Create bonds (use addBond directly since we control placement)
  graph.addBond(attachToId, cId, 1);
  graph.addBond(cId, oId, 2); // Double bond

  return { atomIds: [cId, oId], entryAtomId: cId };
}

/**
 * Apply a carboxyl group (-COOH) to an atom
 */
export function applyCarboxyl(graph: MoleculeGraph, attachToId: string): TemplateInstance | null {
  const attachAtom = graph.atoms.get(attachToId);
  if (!attachAtom) return null;

  // Place C relative to attachment point
  const cPosition = placeAtomRelative(attachAtom, "C");
  const cId = graph.addAtom({ element: "C", position: cPosition });
  const cAtom = graph.atoms.get(cId)!;

  // Place O (double bond)
  const o1Position = placeAtomRelative(cAtom, "O", [1, 0, 0]);
  const o1Id = graph.addAtom({ element: "O", position: o1Position });

  // Place O-H (single bond)
  const o2Position = placeAtomRelative(cAtom, "O", [-1, 0, 0]);
  const o2Id = graph.addAtom({ element: "O", position: o2Position });
  const o2Atom = graph.atoms.get(o2Id)!;

  const hPosition = placeAtomRelative(o2Atom, "H");
  const hId = graph.addAtom({ element: "H", position: hPosition });

  // Create bonds (use addBond directly since we control placement)
  graph.addBond(attachToId, cId, 1);
  graph.addBond(cId, o1Id, 2); // Double bond
  graph.addBond(cId, o2Id, 1);
  graph.addBond(o2Id, hId, 1);

  return { atomIds: [cId, o1Id, o2Id, hId], entryAtomId: cId };
}

/**
 * Apply an amine group (-NH2) to an atom
 */
export function applyAmine(graph: MoleculeGraph, attachToId: string): TemplateInstance | null {
  const attachAtom = graph.atoms.get(attachToId);
  if (!attachAtom) return null;

  // Place N relative to attachment point
  const nPosition = placeAtomRelative(attachAtom, "N");
  const nId = graph.addAtom({ element: "N", position: nPosition });
  const nAtom = graph.atoms.get(nId)!;

  // Place two H atoms
  const h1Position = placeAtomRelative(nAtom, "H", [1, 0, 0]);
  const h1Id = graph.addAtom({ element: "H", position: h1Position });

  const h2Position = placeAtomRelative(nAtom, "H", [-1, 0, 0]);
  const h2Id = graph.addAtom({ element: "H", position: h2Position });

  // Create bonds (use addBond directly since we control placement)
  graph.addBond(attachToId, nId, 1);
  graph.addBond(nId, h1Id, 1);
  graph.addBond(nId, h2Id, 1);

  return { atomIds: [nId, h1Id, h2Id], entryAtomId: nId };
}

/**
 * Apply an amide group (-CONH2) to an atom
 */
export function applyAmide(graph: MoleculeGraph, attachToId: string): TemplateInstance | null {
  const attachAtom = graph.atoms.get(attachToId);
  if (!attachAtom) return null;

  // Place C relative to attachment point
  const cPosition = placeAtomRelative(attachAtom, "C");
  const cId = graph.addAtom({ element: "C", position: cPosition });
  const cAtom = graph.atoms.get(cId)!;

  // Place O (double bond)
  const oPosition = placeAtomRelative(cAtom, "O", [1, 0, 0]);
  const oId = graph.addAtom({ element: "O", position: oPosition });

  // Place N
  const nPosition = placeAtomRelative(cAtom, "N", [-1, 0, 0]);
  const nId = graph.addAtom({ element: "N", position: nPosition });
  const nAtom = graph.atoms.get(nId)!;

  // Place two H atoms on N
  const h1Position = placeAtomRelative(nAtom, "H", [0, 1, 0]);
  const h1Id = graph.addAtom({ element: "H", position: h1Position });

  const h2Position = placeAtomRelative(nAtom, "H", [0, -1, 0]);
  const h2Id = graph.addAtom({ element: "H", position: h2Position });

  // Create bonds (use addBond directly since we control placement)
  graph.addBond(attachToId, cId, 1);
  graph.addBond(cId, oId, 2); // Double bond
  graph.addBond(cId, nId, 1);
  graph.addBond(nId, h1Id, 1);
  graph.addBond(nId, h2Id, 1);

  return { atomIds: [cId, oId, nId, h1Id, h2Id], entryAtomId: cId };
}

