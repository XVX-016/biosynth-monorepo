import { describe, it, expect, beforeEach } from "vitest";
import { MoleculeGraph } from "@biosynth/engine";
import { applyHydroxyl, applyMethyl, applyCarbonyl, applyAmine } from "../groups";
import { useMoleculeStore } from "../../store/moleculeStore";

describe("groups", () => {
  let graph: MoleculeGraph;

  beforeEach(() => {
    useMoleculeStore.getState().reset();
    graph = new MoleculeGraph();
  });

  it("applies hydroxyl group", () => {
    const cId = graph.addAtom({ element: "C", position: [0, 0, 0] });
    const result = applyHydroxyl(graph, cId);

    expect(result).not.toBeNull();
    expect(result!.atomIds.length).toBe(2); // O and H
    expect(result!.entryAtomId).toBeDefined();
    expect(graph.bonds.size).toBe(2); // C-O and O-H
  });

  it("applies methyl group", () => {
    const cId = graph.addAtom({ element: "C", position: [0, 0, 0] });
    const result = applyMethyl(graph, cId);

    expect(result).not.toBeNull();
    expect(result!.atomIds.length).toBe(4); // C and 3 H
    expect(graph.bonds.size).toBe(4); // C-C and 3 C-H
  });

  it("applies carbonyl group", () => {
    const cId = graph.addAtom({ element: "C", position: [0, 0, 0] });
    const result = applyCarbonyl(graph, cId);

    expect(result).not.toBeNull();
    expect(result!.atomIds.length).toBe(2); // C and O
    expect(graph.bonds.size).toBe(2); // C-C and C=O
  });

  it("applies amine group", () => {
    const cId = graph.addAtom({ element: "C", position: [0, 0, 0] });
    const result = applyAmine(graph, cId);

    expect(result).not.toBeNull();
    expect(result!.atomIds.length).toBe(3); // N and 2 H
    expect(graph.bonds.size).toBe(3); // C-N and 2 N-H
  });

  it("returns null for invalid attachment atom", () => {
    const result = applyHydroxyl(graph, "invalid-id");
    expect(result).toBeNull();
  });
});

