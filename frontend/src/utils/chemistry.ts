import type { Element } from '@biosynth/engine'

export const ELEMENT_RADII: Record<string, number> = {
  H: 0.55,
  C: 1,
  O: 0.9,
  N: 0.95,
  F: 0.9,
  S: 1.2,
  P: 1.15,
  Cl: 1.1,
  Br: 1.25,
  I: 1.4,
}

export const ELEMENT_BASE_COLORS: Record<string, number> = {
  H: 0xf6f7f8,
  C: 0x9da3ae,
  O: 0x6b8fa3,
  N: 0x8b95a1,
  F: 0xa8b5b8,
  S: 0xc4b5a0,
  P: 0xa8a39d,
  Cl: 0x9fa8b3,
  Br: 0x8b8b8b,
  I: 0x7a7a7a,
}

export const ELEMENT_ACCENT_COLORS: Record<string, number> = {
  H: 0xffffff,
  C: 0x2b2e33,
  O: 0x4a6fa5,
  N: 0x5a6b7a,
  F: 0x7a8b8f,
  S: 0xb5a085,
  P: 0x8b7a6b,
  Cl: 0x6b7a8b,
  Br: 0x6b6b6b,
  I: 0x5a5a5a,
}

export const VALENCE_MAP: Partial<Record<Element, number>> = {
  H: 1,
  C: 4,
  N: 3,
  O: 2,
  F: 1,
  P: 5,
  S: 6,
  Cl: 1,
  Br: 1,
  I: 1,
}

export const ELECTRONEGATIVITY: Partial<Record<Element, number>> = {
  H: 2.2,
  C: 2.55,
  N: 3.04,
  O: 3.44,
  F: 3.98,
  P: 2.19,
  S: 2.58,
  Cl: 3.16,
  Br: 2.96,
  I: 2.66,
}

export function clamp(value: number, min: number, max: number) {
  return Math.min(max, Math.max(min, value))
}

export function interpolateColor(start: number, end: number, t: number): number {
  const sr = (start >> 16) & 0xff
  const sg = (start >> 8) & 0xff
  const sb = start & 0xff
  const er = (end >> 16) & 0xff
  const eg = (end >> 8) & 0xff
  const eb = end & 0xff

  const rr = Math.round(sr + (er - sr) * t)
  const rg = Math.round(sg + (eg - sg) * t)
  const rb = Math.round(sb + (eb - sb) * t)

  return (rr << 16) | (rg << 8) | rb
}

export const ATOMIC_MASSES: Record<string, number> = {
  H: 1.008,
  C: 12.011,
  N: 14.007,
  O: 15.999,
  F: 18.998,
  P: 30.974,
  S: 32.06,
  Cl: 35.45,
  Br: 79.904,
  I: 126.9,
}

export const BOND_LENGTHS: Record<string, number> = {
  'C-C:1': 1.54,
  'C-C:2': 1.34,
  'C-C:3': 1.20,
  'C-H:1': 1.09,
  'C-O:1': 1.43,
  'C-O:2': 1.23,
  'C-N:1': 1.47,
  'C-N:2': 1.30,
  'C-F:1': 1.35,
  'C-Cl:1': 1.77,
  'C-Br:1': 1.94,
  'C-I:1': 2.14,
  'N-H:1': 1.01,
  'O-H:1': 0.96,
  'S-H:1': 1.34,
}

export const IR_FREQUENCY_WINDOWS: Record<
  string,
  { min: number; max: number; label: string }
> = {
  'C=O': { min: 1680, max: 1750, label: 'Carbonyl stretch' },
  'O-H': { min: 3200, max: 3600, label: 'Hydroxyl stretch' },
  'C-H': { min: 2850, max: 3000, label: 'sp3 C-H stretch' },
  'N-H': { min: 3300, max: 3500, label: 'Amine stretch' },
  'C≡C': { min: 2100, max: 2260, label: 'Alkyne stretch' },
  'C≡N': { min: 2240, max: 2260, label: 'Nitrile stretch' },
  'C=C': { min: 1600, max: 1680, label: 'Alkene stretch' },
  'Ar': { min: 1450, max: 1600, label: 'Aromatic ring stretch' },
}

export const NMR_SHIFT_REFERENCES: Record<
  string,
  { min: number; max: number; multiplicity?: string }
> = {
  'alkyl': { min: 0.8, max: 1.8 },
  'allylic': { min: 2.0, max: 3.0 },
  'aromatic': { min: 6.5, max: 8.5 },
  'aldehyde': { min: 9.0, max: 10.5 },
  'carboxylic': { min: 10.0, max: 13.0 },
  'alcohol': { min: 3.0, max: 4.5 },
}

export function bondKey(a: string, b: string, order: number): string {
  const sorted = [a, b].sort().join('-')
  return `${sorted}:${order}`
}

export function approximateBondLength(
  elA: string,
  elB: string,
  order: number
): number {
  const key = bondKey(elA, elB, order)
  return BOND_LENGTHS[key] ?? 1.5
}

export function getAtomicMass(element: string): number {
  return ATOMIC_MASSES[element] ?? 12
}



