import { describe, it, expect, beforeEach } from "vitest";
import { loadTemplate, placeTemplate, getTemplates, attachTemplateToAtom } from "../templateLoader";
import { useMoleculeStore } from "../../store/moleculeStore";
import { MoleculeGraph } from "@biosynth/engine";

describe("templateLoader", () => {
  beforeEach(() => {
    useMoleculeStore.getState().reset();
  });

  it("loads all templates", () => {
    const templates = getTemplates();
    expect(templates.length).toBeGreaterThan(0);
    expect(templates.length).toBe(4);
  });

  it("loads a specific template", () => {
    const t = loadTemplate("water");
    expect(t).not.toBeNull();
    expect(t!.atoms.length).toBe(3);
    expect(t!.bonds.length).toBe(2);
  });

  it("returns null for invalid template id", () => {
    const t = loadTemplate("invalid");
    expect(t).toBeNull();
  });

  it("places a template into the graph", () => {
    const t = loadTemplate("methane")!;
    const created = placeTemplate(t);

    const molecule = useMoleculeStore.getState().currentMolecule;
    expect(molecule).not.toBeNull();
    expect(created.length).toBe(5); // CH4 = 1C + 4H
    expect(molecule!.atoms.size).toBe(5);
    expect(molecule!.bonds.size).toBe(4);
  });

  it("places template with offset", () => {
    const t = loadTemplate("water")!;
    const created = placeTemplate(t, { x: 5, y: 5, z: 5 });

    const molecule = useMoleculeStore.getState().currentMolecule;
    expect(molecule).not.toBeNull();
    
    // Check that first atom is at offset position
    const firstAtomId = created[0];
    const firstAtom = molecule!.atoms.get(firstAtomId);
    expect(firstAtom).toBeDefined();
    expect(firstAtom!.position[0]).toBeCloseTo(5, 1);
    expect(firstAtom!.position[1]).toBeCloseTo(5, 1);
    expect(firstAtom!.position[2]).toBeCloseTo(5, 1);
  });

  it("appends to existing molecule", () => {
    // Create initial molecule with one atom
    const initialMolecule = new MoleculeGraph();
    initialMolecule.addAtom({ element: "C", position: [10, 10, 10] });
    useMoleculeStore.getState().setMolecule(initialMolecule);

    // Place template
    const t = loadTemplate("water")!;
    const created = placeTemplate(t);

    const molecule = useMoleculeStore.getState().currentMolecule;
    expect(molecule!.atoms.size).toBe(4); // 1 initial + 3 from water
    expect(created.length).toBe(3);
  });

  it("creates molecule if none exists", () => {
    expect(useMoleculeStore.getState().currentMolecule).toBeNull();
    
    const t = loadTemplate("benzene")!;
    placeTemplate(t);

    const molecule = useMoleculeStore.getState().currentMolecule;
    expect(molecule).not.toBeNull();
    expect(molecule!.atoms.size).toBe(6);
    expect(molecule!.bonds.size).toBe(6);
  });

  it("attaches template to atom", () => {
    // Create a molecule with one atom
    const molecule = new MoleculeGraph();
    const atomId = molecule.addAtom({ element: "C", position: [0, 0, 0] });
    useMoleculeStore.getState().setMolecule(molecule);

    const result = attachTemplateToAtom("water", atomId);
    expect(result).not.toBeNull();
    expect(result!.length).toBe(3); // Water has 3 atoms

    const updatedMolecule = useMoleculeStore.getState().currentMolecule;
    expect(updatedMolecule!.atoms.size).toBe(4); // 1 C + 3 from water
    expect(updatedMolecule!.bonds.size).toBeGreaterThan(2); // Original bonds + attachment bond
  });

  it("returns null for invalid template id", () => {
    const molecule = new MoleculeGraph();
    const atomId = molecule.addAtom({ element: "C", position: [0, 0, 0] });
    useMoleculeStore.getState().setMolecule(molecule);

    const result = attachTemplateToAtom("invalid", atomId);
    expect(result).toBeNull();
  });

  it("returns null for invalid atom id", () => {
    const molecule = new MoleculeGraph();
    useMoleculeStore.getState().setMolecule(molecule);

    const result = attachTemplateToAtom("water", "invalid-id");
    expect(result).toBeNull();
  });
});

