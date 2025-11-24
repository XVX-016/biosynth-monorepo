import type { Atom, Bond, Molecule } from '../types/molecule'

// Simple covalent radii (Å) and valences (typical)
const covalentRadius: Record<string, number> = {
  H: 0.31, C: 0.76, N: 0.71, O: 0.66, F: 0.57,
  P: 1.07, S: 1.05, Cl: 1.02, Br: 1.2, I: 1.39
}

const typicalValence: Record<string, number[]> = {
  H: [1],
  C: [4],
  N: [3, 4],
  O: [2],
  S: [2, 4, 6],
  P: [3, 5],
  F: [1],
  Cl: [1],
  Br: [1],
  I: [1],
}

export function maxAllowedBonds(atom: Atom): number {
  const vals = typicalValence[atom.element]
  if (!vals) return 4
  return vals[0]
}

export function distance(a: Atom, b: Atom): number {
  const dx = a.position[0] - b.position[0]
  const dy = a.position[1] - b.position[1]
  const dz = a.position[2] - b.position[2]
  return Math.sqrt(dx*dx + dy*dy + dz*dz)
}

export function shouldAutoBond(a: Atom, b: Atom, mol?: Molecule, tolerance = 0.45): boolean {
  const rA = covalentRadius[a.element] || 0.8
  const rB = covalentRadius[b.element] || 0.8
  const threshold = rA + rB + tolerance
  const d = distance(a, b)
  if (d > threshold) return false

  // check current valence usage
  const molBonds = mol?.bonds || []
  const countA = molBonds.filter(x => x.atom1 === a.id || x.atom2 === a.id).length
  const countB = molBonds.filter(x => x.atom1 === b.id || x.atom2 === b.id).length
  if (countA >= maxAllowedBonds(a) || countB >= maxAllowedBonds(b)) return false

  return true
}

// Apply autoBond across a whole molecule: returns new bonds to add
export function computeAutoBonds(mol: Molecule, tolerance = 0.45): Array<{a: string, b: string}> {
  const candidates: Array<{a: string, b: string}> = []
  const existingBondPairs = new Set<string>()
  
  // Track existing bonds
  mol.bonds.forEach(bond => {
    const key1 = `${bond.atom1}-${bond.atom2}`
    const key2 = `${bond.atom2}-${bond.atom1}`
    existingBondPairs.add(key1)
    existingBondPairs.add(key2)
  })
  
  for (let i = 0; i < mol.atoms.length; i++) {
    for (let j = i + 1; j < mol.atoms.length; j++) {
      const a = mol.atoms[i]
      const b = mol.atoms[j]
      
      // Skip if bond already exists
      const bondKey = `${a.id}-${b.id}`
      if (existingBondPairs.has(bondKey)) continue
      
      if (shouldAutoBond(a, b, mol, tolerance)) {
        candidates.push({ a: a.id, b: b.id })
      }
    }
  }
  return candidates
}

