import { describe, it, expect } from "vitest";
import { MoleculeGraph } from "../src/MoleculeGraph";
import { Atom } from "../src/types";
import { Bond } from "../src/types";

describe("MoleculeGraph", () => {
  it("adds atoms and bonds and neighbors", () => {
    const g = new MoleculeGraph();
    const a1: Atom = { id: "a1", element: "C", position: [0, 0, 0] };
    const a2: Atom = { id: "a2", element: "H", position: [1, 0, 0] };
    g.addAtom({ element: "C", position: [0, 0, 0] });
    g.addAtom({ element: "H", position: [1, 0, 0] });
    
    const atomIds = Array.from(g.atoms.keys());
    const bondId = g.addBond(atomIds[0], atomIds[1], 1);
    expect(bondId).toBeTruthy();
    
    const n = g.getNeighbors(atomIds[0]);
    expect(n).toContain(atomIds[1]);
  });

  it("computes degree correctly", () => {
    const g = new MoleculeGraph();
    const cId = g.addAtom({ element: "C", position: [0, 0, 0] });
    const h1Id = g.addAtom({ element: "H", position: [1, 0, 0] });
    const h2Id = g.addAtom({ element: "H", position: [-1, 0, 0] });
    
    g.addBond(cId, h1Id, 1);
    g.addBond(cId, h2Id, 1);
    
    expect(g.computeDegree(cId)).toBe(2);
    expect(g.computeDegree(h1Id)).toBe(1);
  });

  it("calculates formula correctly", () => {
    const g = new MoleculeGraph();
    const cId = g.addAtom({ element: "C", position: [0, 0, 0] });
    const h1Id = g.addAtom({ element: "H", position: [1, 0, 0] });
    const h2Id = g.addAtom({ element: "H", position: [-1, 0, 0] });
    const h3Id = g.addAtom({ element: "H", position: [0, 1, 0] });
    const h4Id = g.addAtom({ element: "H", position: [0, -1, 0] });
    
    g.addBond(cId, h1Id, 1);
    g.addBond(cId, h2Id, 1);
    g.addBond(cId, h3Id, 1);
    g.addBond(cId, h4Id, 1);
    
    const formula = g.getFormula();
    expect(formula).toContain("C");
    expect(formula).toContain("H");
  });
});

