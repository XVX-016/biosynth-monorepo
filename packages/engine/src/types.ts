// Core types for the molecular engine

export type Element = "H" | "C" | "N" | "O" | "F" | "S" | "P" | "Cl" | "Br" | "I";

export interface Atom {
  id: string;
  element: Element;
  position: [number, number, number];
  charge?: number;
}

export interface Bond {
  id: string;
  a1: string; // atom ID
  a2: string; // atom ID
  order: number; // 1 (single), 2 (double), 3 (triple)
}

export interface MoleculeState {
  atoms: Map<string, Atom>;
  bonds: Map<string, Bond>;
}

// Element properties
export const ELEMENT_COLORS: Record<Element, string> = {
  H: "#3b82f6", // Blue
  C: "#1e293b", // Dark gray/black
  N: "#8b5cf6", // Purple
  O: "#ef4444", // Red
  F: "#10b981", // Green
  S: "#f59e0b", // Orange
  P: "#f97316", // Orange-red
  Cl: "#84cc16", // Lime
  Br: "#991b1b", // Dark red
  I: "#7c3aed", // Violet
};

export const ELEMENT_RADII: Record<Element, number> = {
  H: 0.31,
  C: 0.77,
  N: 0.75,
  O: 0.73,
  F: 0.72,
  S: 1.02,
  P: 1.06,
  Cl: 0.99,
  Br: 1.14,
  I: 1.33,
};

export const COVALENT_BOND_LENGTHS: Record<string, number> = {
  "C-C": 1.54,
  "C-H": 1.09,
  "C-O": 1.43,
  "C-N": 1.47,
  "O-H": 0.96,
  "N-H": 1.01,
  "C-F": 1.35,
  "C-Cl": 1.77,
  "C-Br": 1.94,
  "C-I": 2.14,
  "O-O": 1.48,
  "N-N": 1.45,
};

export function getBondLength(e1: Element, e2: Element): number {
  const key1 = `${e1}-${e2}`;
  const key2 = `${e2}-${e1}`;
  return COVALENT_BOND_LENGTHS[key1] || COVALENT_BOND_LENGTHS[key2] || 1.5;
}

