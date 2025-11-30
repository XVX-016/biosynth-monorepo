import type { MoleculeGraph } from '@biosynth/engine'

export interface StructureIssueSummary {
  unbondedAtoms: number
  hasInstabilityRisk: boolean
  warnings: string[]
}

/**
 * Very lightweight structural checks so the UI can surface actionable issues.
 * This intentionally runs on every render, so keep it O(n).
 */
export function analyzeStructure(molecule: MoleculeGraph | null): StructureIssueSummary {
  if (!molecule) {
    return {
      unbondedAtoms: 0,
      hasInstabilityRisk: false,
      warnings: ['No molecule loaded'],
    }
  }

  let unbondedAtoms = 0
  const bondedAtoms = new Set<string>()

  molecule.bonds.forEach((bond) => {
    bondedAtoms.add(bond.a1)
    bondedAtoms.add(bond.a2)
  })

  molecule.atoms.forEach((atom) => {
    if (!bondedAtoms.has(atom.id)) {
      unbondedAtoms += 1
    }
  })

  const warnings: string[] = []
  if (unbondedAtoms > 0) {
    warnings.push(`${unbondedAtoms} atom${unbondedAtoms === 1 ? '' : 's'} have no bonds`)
  }

  if (molecule.bonds.size === 0 && molecule.atoms.size > 3) {
    warnings.push('Structure may be unstable without bonds')
  }

  return {
    unbondedAtoms,
    hasInstabilityRisk: warnings.length > 0,
    warnings: warnings.length ? warnings : ['Structure looks stable'],
  }
}


