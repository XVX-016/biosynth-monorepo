// src/templates/atoms.ts
import { Element } from "@biosynth/engine";

/**
 * Atomic templates - preconfigured atoms with default properties
 * These can be used to quickly add atoms to molecules
 */
export const atomTemplates: Record<string, { element: Element; position: [number, number, number]; charge?: number }> = {
  C: { element: "C", position: [0, 0, 0], charge: 0 },
  H: { element: "H", position: [0, 0, 0], charge: 0 },
  O: { element: "O", position: [0, 0, 0], charge: 0 },
  N: { element: "N", position: [0, 0, 0], charge: 0 },
  S: { element: "S", position: [0, 0, 0], charge: 0 },
  P: { element: "P", position: [0, 0, 0], charge: 0 },
  F: { element: "F", position: [0, 0, 0], charge: 0 },
  Cl: { element: "Cl", position: [0, 0, 0], charge: 0 },
  Br: { element: "Br", position: [0, 0, 0], charge: 0 },
  I: { element: "I", position: [0, 0, 0], charge: 0 },
};

export type AtomTemplateName = keyof typeof atomTemplates;

