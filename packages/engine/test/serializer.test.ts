import { describe, it, expect } from "vitest";
import { MoleculeGraph } from "../src/MoleculeGraph";
import { MoleculeSerializer } from "../src/MoleculeSerializer";

describe("MoleculeSerializer", () => {
  it("serializes and deserializes JSON correctly", () => {
    const graph = new MoleculeGraph();
    const cId = graph.addAtom({ element: "C", position: [0, 0, 0] });
    const hId = graph.addAtom({ element: "H", position: [1, 0, 0] });
    graph.addBond(cId, hId, 1);

    const json = MoleculeSerializer.toJSON(graph);
    expect(json.atoms).toHaveLength(2);
    expect(json.bonds).toHaveLength(1);

    const restored = MoleculeSerializer.fromJSON(json);
    expect(restored.atoms.size).toBe(2);
    expect(restored.bonds.size).toBe(1);
  });

  it("converts to SMILES (placeholder)", () => {
    const graph = new MoleculeGraph();
    graph.addAtom({ element: "C", position: [0, 0, 0] });

    const smiles = MoleculeSerializer.toSMILES(graph);
    expect(smiles).toBe("C");
  });

  it("returns empty string for empty molecule", () => {
    const graph = new MoleculeGraph();
    const smiles = MoleculeSerializer.toSMILES(graph);
    expect(smiles).toBe("");
  });

  it("fromSMILES returns null (placeholder)", () => {
    const result = MoleculeSerializer.fromSMILES("C");
    expect(result).toBeNull();
  });
});

