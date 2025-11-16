import { describe, it, expect } from "vitest";
import { createBenzene, createMethane, createWater, createEthanol } from "../molecules";

describe("molecules", () => {
  it("creates benzene ring", () => {
    const molecule = createBenzene();
    expect(molecule.atoms.size).toBe(12); // 6 C + 6 H
    expect(molecule.bonds.size).toBe(12); // 6 ring bonds + 6 C-H bonds
  });

  it("creates methane", () => {
    const molecule = createMethane();
    expect(molecule.atoms.size).toBe(5); // 1 C + 4 H
    expect(molecule.bonds.size).toBe(4); // 4 C-H bonds
  });

  it("creates water", () => {
    const molecule = createWater();
    expect(molecule.atoms.size).toBe(3); // 1 O + 2 H
    expect(molecule.bonds.size).toBe(2); // 2 O-H bonds
  });

  it("creates ethanol", () => {
    const molecule = createEthanol();
    expect(molecule.atoms.size).toBe(9); // 2 C + 1 O + 6 H
    expect(molecule.bonds.size).toBeGreaterThan(5);
  });

  it("benzene has correct formula", () => {
    const molecule = createBenzene();
    const formula = molecule.getFormula();
    expect(formula).toContain("C");
    expect(formula).toContain("H");
  });
});

