import { describe, it, expect } from 'vitest';
import { getValence, shouldFormBond, suggestBondOrder } from '../valence';
import { Atom } from '@biosynth/engine';

describe('valence helpers', () => {
  it('returns known valences', () => {
    expect(getValence('C')).toBe(4);
    expect(getValence('O')).toBe(2);
    expect(getValence('H')).toBe(1);
    expect(getValence('N')).toBe(3);
    expect(getValence('Xx')).toBe(0);
  });

  it('suggest bond order based on remaining valence', () => {
    const a: Atom = {
      id: 'a1',
      element: 'C',
      position: [0, 0, 0],
    };
    const b: Atom = {
      id: 'a2',
      element: 'C',
      position: [1.4, 0, 0],
    };
    // Both have full valence available - defaults to single bond (conservative)
    expect(suggestBondOrder(a, b, [])).toBe(1);
    
    // When both atoms have limited remaining valence, suggest higher order
    const bonds = [
      { a1: 'a1', a2: 'other1', order: 2 },
      { a1: 'a2', a2: 'other2', order: 2 },
    ];
    // Both have 2 remaining, so can form double bond
    expect(suggestBondOrder(a, b, bonds)).toBe(2);
    
    // Triple bond case: both have exactly 3 remaining
    const tripleBonds = [
      { a1: 'a1', a2: 'other1', order: 1 },
      { a1: 'a2', a2: 'other2', order: 1 },
    ];
    expect(suggestBondOrder(a, b, tripleBonds)).toBe(3);
  });

  it('should form bond when in threshold and valence available', () => {
    const a: Atom = {
      id: 'a1',
      element: 'C',
      position: [0, 0, 0],
    };
    const b: Atom = {
      id: 'a2',
      element: 'H',
      position: [0.9, 0, 0],
    };
    expect(shouldFormBond(a, b, [], 1.2)).toBe(true);
  });

  it('should not form bond when out of threshold', () => {
    const a: Atom = {
      id: 'a1',
      element: 'C',
      position: [0, 0, 0],
    };
    const b: Atom = {
      id: 'a2',
      element: 'C',
      position: [3.0, 0, 0],
    };
    expect(shouldFormBond(a, b, [], 1.6)).toBe(false);
  });

  it('should not form bond when valence exhausted', () => {
    const a: Atom = {
      id: 'a1',
      element: 'H',
      position: [0, 0, 0],
    };
    const b: Atom = {
      id: 'a2',
      element: 'H',
      position: [0.9, 0, 0],
    };
    // Both H atoms already bonded (valence 1)
    const bonds = [
      { a1: 'a1', a2: 'other', order: 1 },
      { a2: 'a2', a1: 'other2', order: 1 },
    ];
    expect(shouldFormBond(a, b, bonds, 1.2)).toBe(false);
  });
});

