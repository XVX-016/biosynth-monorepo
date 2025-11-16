import { describe, it, expect, beforeEach } from 'vitest';
import { MoleculeGraph } from '@biosynth/engine';
import { createBondSafe, removeBondSafe, hasBondBetween, findBondBetween } from '../bonds';

describe('bond kernel', () => {
  let molecule: MoleculeGraph;

  beforeEach(() => {
    molecule = new MoleculeGraph();
    // add two atoms
    molecule.addAtom({ element: 'C', position: [0, 0, 0] });
    molecule.addAtom({ element: 'C', position: [1.4, 0, 0] });
  });

  it('creates a bond when atoms are close and valence allows', () => {
    const atoms = Array.from(molecule.atoms.values());
    const [a, b] = atoms;
    const bondId = createBondSafe(molecule, a.id, b.id);
    expect(bondId).not.toBeNull();
    expect(hasBondBetween(molecule, a.id, b.id)).toBe(true);
    expect(findBondBetween(molecule, a.id, b.id)).toBeDefined();
  });

  it('does not create duplicate bond', () => {
    const atoms = Array.from(molecule.atoms.values());
    const [a, b] = atoms;
    const b1 = createBondSafe(molecule, a.id, b.id);
    expect(b1).not.toBeNull();
    const b2 = createBondSafe(molecule, a.id, b.id);
    expect(b2).toBe(b1); // Should return existing bond ID
  });

  it('removes bond from molecule', () => {
    const atoms = Array.from(molecule.atoms.values());
    const [a, b] = atoms;
    const bondId = createBondSafe(molecule, a.id, b.id);
    expect(bondId).not.toBeNull();
    
    const removed = removeBondSafe(molecule, bondId!);
    expect(removed).toBe(true);
    expect(hasBondBetween(molecule, a.id, b.id)).toBe(false);
  });

  it('does not create bond when atoms are too far apart', () => {
    const atoms = Array.from(molecule.atoms.values());
    const [a] = atoms;
    // Add atom far away
    const farAtomId = molecule.addAtom({ element: 'C', position: [10, 10, 10] });
    const bondId = createBondSafe(molecule, a.id, farAtomId);
    expect(bondId).toBeNull();
  });

  it('finds bond between atoms', () => {
    const atoms = Array.from(molecule.atoms.values());
    const [a, b] = atoms;
    const bondId = createBondSafe(molecule, a.id, b.id);
    expect(bondId).not.toBeNull();
    
    const bond = findBondBetween(molecule, a.id, b.id);
    expect(bond).toBeDefined();
    expect(bond?.id).toBe(bondId);
  });
});

