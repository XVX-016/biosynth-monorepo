import { describe, it, expect } from "vitest";
import { MoleculeGraph } from "../src/MoleculeGraph";

describe("Bond Management", () => {
  it("adds bonds between atoms", () => {
    const graph = new MoleculeGraph();
    const cId = graph.addAtom({ element: "C", position: [0, 0, 0] });
    const hId = graph.addAtom({ element: "H", position: [1, 0, 0] });

    const bondId = graph.addBond(cId, hId, 1);
    expect(bondId).toBeTruthy();
    expect(graph.bonds.size).toBe(1);
  });

  it("prevents duplicate bonds", () => {
    const graph = new MoleculeGraph();
    const cId = graph.addAtom({ element: "C", position: [0, 0, 0] });
    const hId = graph.addAtom({ element: "H", position: [1, 0, 0] });

    const bondId1 = graph.addBond(cId, hId, 1);
    const bondId2 = graph.addBond(cId, hId, 1);
    const bondId3 = graph.addBond(hId, cId, 1); // Reverse order

    expect(bondId1).toBeTruthy();
    expect(bondId2).toBe(bondId1); // Should return existing bond
    expect(bondId3).toBe(bondId1); // Should return existing bond
    expect(graph.bonds.size).toBe(1);
  });

  it("updates graph when bond is added", () => {
    const graph = new MoleculeGraph();
    const cId = graph.addAtom({ element: "C", position: [0, 0, 0] });
    const h1Id = graph.addAtom({ element: "H", position: [1, 0, 0] });
    const h2Id = graph.addAtom({ element: "H", position: [-1, 0, 0] });

    graph.addBond(cId, h1Id, 1);
    graph.addBond(cId, h2Id, 1);

    expect(graph.bonds.size).toBe(2);
    expect(graph.getNeighbors(cId)).toHaveLength(2);
    expect(graph.getNeighbors(h1Id)).toHaveLength(1);
  });

  it("removes bonds when atom is removed", () => {
    const graph = new MoleculeGraph();
    const cId = graph.addAtom({ element: "C", position: [0, 0, 0] });
    const hId = graph.addAtom({ element: "H", position: [1, 0, 0] });

    graph.addBond(cId, hId, 1);
    expect(graph.bonds.size).toBe(1);

    graph.removeAtom(cId);
    expect(graph.bonds.size).toBe(0);
  });

  it("updates neighbors when bond is removed", () => {
    const graph = new MoleculeGraph();
    const cId = graph.addAtom({ element: "C", position: [0, 0, 0] });
    const hId = graph.addAtom({ element: "H", position: [1, 0, 0] });

    const bondId = graph.addBond(cId, hId, 1);
    expect(graph.getNeighbors(cId)).toContain(hId);

    graph.removeBond(bondId!);
    expect(graph.getNeighbors(cId)).not.toContain(hId);
  });
});

