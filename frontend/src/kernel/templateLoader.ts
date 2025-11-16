import { MoleculeGraph } from "@biosynth/engine";
import { useMoleculeStore } from "../store/moleculeStore";
import { pushState } from "../store/historyStore";
import TEMPLATES from "../templates/templates.index";
import { suggestBondVector } from "../templates/placement";
import { getBondLength } from "@biosynth/engine";
import { resolveClashes, planarizeAromaticRings, relaxGeometry } from "./geometryCleanup";

export interface TemplateData {
  atoms: { element: string; coords: { x: number; y: number; z: number } }[];
  bonds: { a: number; b: number; order: number }[];
}

export function getTemplates() {
  return TEMPLATES;
}

export function loadTemplate(id: string): TemplateData | null {
  const t = TEMPLATES.find((t) => t.id === id);
  return t ? t.data : null;
}

/**
 * Places a template into the current molecule.
 * Returns list of created atom IDs for selection & placement logic.
 */
export function placeTemplate(template: TemplateData, offset = { x: 0, y: 0, z: 0 }) {
  const store = useMoleculeStore.getState();
  
  // Get or create molecule
  let molecule = store.currentMolecule;
  if (!molecule) {
    molecule = new MoleculeGraph();
  } else {
    // Clone to avoid mutating the current molecule directly
    molecule = molecule.clone();
  }

  const createdAtomIds: string[] = [];

  // 1. create atoms
  template.atoms.forEach((atom) => {
    const id = molecule.addAtom({
      element: atom.element as any,
      position: [
        atom.coords.x + offset.x,
        atom.coords.y + offset.y,
        atom.coords.z + offset.z
      ]
    });
    createdAtomIds.push(id);
  });

  // 2. create bonds
  template.bonds.forEach((b) => {
    molecule.addBond(createdAtomIds[b.a], createdAtomIds[b.b], b.order);
  });

  // Apply geometry cleanup
  resolveClashes(molecule);
  planarizeAromaticRings(molecule);
  relaxGeometry(molecule);

  // Update store and push to history
  store.setMolecule(molecule);
  pushState();
  
  return createdAtomIds;
}

/**
 * Attach a template to an existing atom in the molecule.
 * Computes direction vector and aligns template's root atom with attachment point.
 */
export function attachTemplateToAtom(templateId: string, atomId: string): string[] | null {
  const store = useMoleculeStore.getState();
  const molecule = store.currentMolecule;
  
  if (!molecule) {
    // Create new molecule if none exists
    const newMolecule = new MoleculeGraph();
    useMoleculeStore.getState().setMolecule(newMolecule);
    return attachTemplateToAtom(templateId, atomId);
  }

  const attachAtom = molecule.atoms.get(atomId);
  if (!attachAtom) return null;

  const template = loadTemplate(templateId);
  if (!template || template.atoms.length === 0) return null;

  // Clone molecule to avoid direct mutation
  const cloned = molecule.clone();

  // Get direction vector for placement
  const direction = suggestBondVector(cloned, atomId);
  const bondLength = getBondLength(attachAtom.element, template.atoms[0].element as any);

  // Root atom of template (first atom)
  const rootAtomTemplate = template.atoms[0];
  
  // Calculate offset to align root atom with attachment point
  // Root atom should be placed at: attachAtom.position + direction * bondLength
  const rootPosition: [number, number, number] = [
    attachAtom.position[0] + direction.x * bondLength,
    attachAtom.position[1] + direction.y * bondLength,
    attachAtom.position[2] + direction.z * bondLength,
  ];

  // Calculate template's root atom offset from origin
  const rootOffset = {
    x: rootAtomTemplate.coords.x,
    y: rootAtomTemplate.coords.y,
    z: rootAtomTemplate.coords.z,
  };

  // Create atoms with proper offset
  const createdAtomIds: string[] = [];
  template.atoms.forEach((atom) => {
    const id = cloned.addAtom({
      element: atom.element as any,
      position: [
        rootPosition[0] + (atom.coords.x - rootOffset.x),
        rootPosition[1] + (atom.coords.y - rootOffset.y),
        rootPosition[2] + (atom.coords.z - rootOffset.z),
      ]
    });
    createdAtomIds.push(id);
  });

  // Create bonds within template
  template.bonds.forEach((b) => {
    cloned.addBond(createdAtomIds[b.a], createdAtomIds[b.b], b.order);
  });

  // Create bond between attachment atom and template root
  const rootAtomId = createdAtomIds[0];
  cloned.addBond(atomId, rootAtomId, 1);

  // Apply geometry cleanup
  resolveClashes(cloned);
  planarizeAromaticRings(cloned);
  relaxGeometry(cloned);

  // Update store and push to history
  store.setMolecule(cloned);
  pushState();

  return createdAtomIds;
}

