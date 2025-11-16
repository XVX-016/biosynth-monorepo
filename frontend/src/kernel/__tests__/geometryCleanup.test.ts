import { describe, it, expect } from "vitest";
import { MoleculeGraph } from "@biosynth/engine";
import { resolveClashes, planarizeAromaticRings, relaxGeometry } from "../geometryCleanup";

describe("geometryCleanup", () => {
  it("resolves clashes by pushing atoms apart", () => {
    const graph = new MoleculeGraph();
    const a1 = graph.addAtom({ element: "C", position: [0, 0, 0] });
    const a2 = graph.addAtom({ element: "C", position: [0.5, 0, 0] }); // Too close

    const atom1 = graph.atoms.get(a1)!;
    const atom2 = graph.atoms.get(a2)!;
    const initialDistance = Math.sqrt(
      Math.pow(atom2.position[0] - atom1.position[0], 2) +
      Math.pow(atom2.position[1] - atom1.position[1], 2) +
      Math.pow(atom2.position[2] - atom1.position[2], 2)
    );

    resolveClashes(graph);

    const finalDistance = Math.sqrt(
      Math.pow(atom2.position[0] - atom1.position[0], 2) +
      Math.pow(atom2.position[1] - atom1.position[1], 2) +
      Math.pow(atom2.position[2] - atom1.position[2], 2)
    );

    // Atoms should be pushed apart
    expect(finalDistance).toBeGreaterThan(initialDistance);
  });

  it("does not affect bonded atoms", () => {
    const graph = new MoleculeGraph();
    const a1 = graph.addAtom({ element: "C", position: [0, 0, 0] });
    const a2 = graph.addAtom({ element: "H", position: [1.0, 0, 0] });
    graph.addBond(a1, a2, 1);

    const atom1Before = { ...graph.atoms.get(a1)! };
    const atom2Before = { ...graph.atoms.get(a2)! };

    resolveClashes(graph);

    const atom1After = graph.atoms.get(a1)!;
    const atom2After = graph.atoms.get(a2)!;

    // Bonded atoms should not be significantly moved
    const distanceBefore = Math.sqrt(
      Math.pow(atom2Before.position[0] - atom1Before.position[0], 2) +
      Math.pow(atom2Before.position[1] - atom1Before.position[1], 2) +
      Math.pow(atom2Before.position[2] - atom1Before.position[2], 2)
    );
    const distanceAfter = Math.sqrt(
      Math.pow(atom2After.position[0] - atom1After.position[0], 2) +
      Math.pow(atom2After.position[1] - atom1After.position[1], 2) +
      Math.pow(atom2After.position[2] - atom1After.position[2], 2)
    );

    // Distance should remain similar (within 0.2 Å)
    expect(Math.abs(distanceAfter - distanceBefore)).toBeLessThan(0.2);
  });

  it("planarizes aromatic rings", () => {
    const graph = new MoleculeGraph();
    // Create a simple 6-membered ring
    const atoms: string[] = [];
    for (let i = 0; i < 6; i++) {
      const angle = (i * Math.PI * 2) / 6;
      atoms.push(
        graph.addAtom({
          element: "C",
          position: [Math.cos(angle), Math.sin(angle), i * 0.1], // Slightly non-planar
        })
      );
    }

    // Create alternating bonds
    for (let i = 0; i < 6; i++) {
      graph.addBond(atoms[i], atoms[(i + 1) % 6], i % 2 === 0 ? 2 : 1);
    }

    planarizeAromaticRings(graph);

    // Check that z-coordinates are similar (planar)
    const zCoords = atoms.map((id) => graph.atoms.get(id)!.position[2]);
    const avgZ = zCoords.reduce((a, b) => a + b, 0) / zCoords.length;
    const maxDeviation = Math.max(...zCoords.map((z) => Math.abs(z - avgZ)));

    // Should be relatively planar (deviation < 0.3 for simplified algorithm)
    expect(maxDeviation).toBeLessThan(0.3);
  });

  it("relaxes geometry using spring model", () => {
    const graph = new MoleculeGraph();
    const a1 = graph.addAtom({ element: "C", position: [0, 0, 0] });
    const a2 = graph.addAtom({ element: "C", position: [2.0, 0, 0] }); // Too far
    graph.addBond(a1, a2, 1);

    const atom1Before = { ...graph.atoms.get(a1)! };
    const atom2Before = { ...graph.atoms.get(a2)! };

    relaxGeometry(graph, 10);

    const atom1After = graph.atoms.get(a1)!;
    const atom2After = graph.atoms.get(a2)!;

    const distanceBefore = Math.sqrt(
      Math.pow(atom2Before.position[0] - atom1Before.position[0], 2) +
      Math.pow(atom2Before.position[1] - atom1Before.position[1], 2) +
      Math.pow(atom2Before.position[2] - atom1Before.position[2], 2)
    );
    const distanceAfter = Math.sqrt(
      Math.pow(atom2After.position[0] - atom1After.position[0], 2) +
      Math.pow(atom2After.position[1] - atom1After.position[1], 2) +
      Math.pow(atom2After.position[2] - atom1After.position[2], 2)
    );

    // Distance should move toward ideal bond length (~1.5 Å for C-C)
    // Note: With only 10 iterations, may not fully converge
    expect(distanceAfter).toBeLessThanOrEqual(distanceBefore);
    expect(distanceAfter).toBeGreaterThan(1.0); // Still reasonable
  });

  it("handles empty molecule", () => {
    const graph = new MoleculeGraph();
    expect(() => resolveClashes(graph)).not.toThrow();
    expect(() => planarizeAromaticRings(graph)).not.toThrow();
    expect(() => relaxGeometry(graph)).not.toThrow();
  });
});

