import type { Molecule } from '../types/molecule'

export function serializeMolecule(mol: Molecule): string {
  return JSON.stringify(mol, null, 2)
}

export function deserializeMolecule(json: string): Molecule {
  try {
    const parsed = JSON.parse(json)
    // basic validation
    parsed.atoms = parsed.atoms || []
    parsed.bonds = parsed.bonds || []
    return parsed as Molecule
  } catch (e) {
    throw new Error('Invalid molecule JSON')
  }
}

