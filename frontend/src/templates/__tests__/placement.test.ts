import { describe, it, expect } from "vitest";
import { MoleculeGraph } from "@biosynth/engine";
import { suggestBondVector, placeTemplate } from "../placement";

describe("placement", () => {
  it("suggests bond vector for atom with no neighbors", () => {
    const graph = new MoleculeGraph();
    const atomId = graph.addAtom({ element: "C", position: [0, 0, 0] });

    const vector = suggestBondVector(graph, atomId);
    expect(vector).toHaveProperty("x");
    expect(vector).toHaveProperty("y");
    expect(vector).toHaveProperty("z");
  });

  it("suggests bond vector away from existing neighbors", () => {
    const graph = new MoleculeGraph();
    const cId = graph.addAtom({ element: "C", position: [0, 0, 0] });
    const hId = graph.addAtom({ element: "H", position: [1, 0, 0] });
    graph.addBond(cId, hId, 1);

    const vector = suggestBondVector(graph, cId);
    // Should point away from the H atom (negative x direction)
    expect(vector.x).toBeLessThan(0);
  });

  it("places template atoms relative to attachment point", () => {
    const graph = new MoleculeGraph();
    const attachId = graph.addAtom({ element: "C", position: [0, 0, 0] });
    const newAtom1Id = graph.addAtom({ element: "H", position: [0, 0, 0] });
    const newAtom2Id = graph.addAtom({ element: "H", position: [0, 0, 0] });

    placeTemplate(graph, [newAtom1Id, newAtom2Id], attachId);

    const newAtom1 = graph.atoms.get(newAtom1Id);
    const newAtom2 = graph.atoms.get(newAtom2Id);
    expect(newAtom1).toBeDefined();
    expect(newAtom2).toBeDefined();
    // Atoms should be moved from origin
    expect(newAtom1!.position[0]).not.toBe(0);
  });

  it("handles empty atom list", () => {
    const graph = new MoleculeGraph();
    const attachId = graph.addAtom({ element: "C", position: [0, 0, 0] });

    expect(() => placeTemplate(graph, [], attachId)).not.toThrow();
  });
});

