const generateId = () => {
  if (typeof crypto !== 'undefined' && typeof crypto.randomUUID === 'function') {
    return crypto.randomUUID()
  }
  return Math.random().toString(36).slice(2)
}
import { MoleculeGraph, ForceField } from '@biosynth/engine'
import type {
  ValidationResult,
  ValidationIssue,
  AutoFix,
} from '../types/validation'
import {
  ELEMENT_RADII,
  approximateBondLength,
  bondKey,
} from '../utils/chemistry'

interface SimpleAtom {
  id: string
  element: string
  position: [number, number, number]
}

interface SimpleBond {
  id: string
  a1: string
  a2: string
  order: number
}

export class MolecularValidator {
  validate(molecule: MoleculeGraph): ValidationResult {
    const atoms = this.toAtoms(molecule)
    const bonds = this.toBonds(molecule)

    const issues: ValidationIssue[] = [
      ...this.checkStericClashes(atoms),
      ...this.checkBondLengths(atoms, bonds),
      ...this.checkAngleStrain(atoms, bonds),
    ]

    const severityPenalty = issues.reduce((penalty, issue) => {
      return penalty + (issue.severity === 'error' ? 10 : 5)
    }, 0)

    const score = Math.max(0, 100 - severityPenalty)

    const suggestions = this.generateSuggestions(issues)
    const autoFixes: AutoFix[] = issues.length
      ? [
          {
            id: 'optimize',
            label: 'Geometry optimization',
            description:
              'Applies a molecular mechanics optimization to relieve clashes and strain.',
            appliesTo: issues.map((issue) => issue.id),
          },
        ]
      : []

    return {
      issues,
      score,
      suggestions,
      autoFixes,
      timestamp: new Date().toISOString(),
    }
  }

  private toAtoms(molecule: MoleculeGraph): SimpleAtom[] {
    const atoms: SimpleAtom[] = []
    molecule.atoms.forEach((atom, id) => {
      atoms.push({
        id,
        element: atom.element,
        position: atom.position,
      })
    })
    return atoms
  }

  private toBonds(molecule: MoleculeGraph): SimpleBond[] {
    const bonds: SimpleBond[] = []
    molecule.bonds.forEach((bond, id) => {
      bonds.push({
        id,
        a1: bond.a1,
        a2: bond.a2,
        order: bond.order,
      })
    })
    return bonds
  }

  private checkStericClashes(atoms: SimpleAtom[]): ValidationIssue[] {
    const issues: ValidationIssue[] = []
    for (let i = 0; i < atoms.length; i++) {
      for (let j = i + 1; j < atoms.length; j++) {
        const atomA = atoms[i]
        const atomB = atoms[j]
        const distance = this.distance(atomA.position, atomB.position)
        const minDistance =
          (ELEMENT_RADII[atomA.element] ?? 1) +
          (ELEMENT_RADII[atomB.element] ?? 1)
        if (distance < minDistance * 0.6) {
          issues.push({
            id: generateId(),
            type: 'steric_clash',
            severity: 'error',
            message: `Steric clash between ${atomA.element}-${atomB.element}`,
            atoms: [atomA.id, atomB.id],
            expectedValue: parseFloat((minDistance * 0.8).toFixed(2)),
            actualValue: parseFloat(distance.toFixed(2)),
            recommendation: 'Increase distance between atoms or optimize geometry.',
          })
        }
      }
    }
    return issues
  }

  private checkBondLengths(
    atoms: SimpleAtom[],
    bonds: SimpleBond[]
  ): ValidationIssue[] {
    const issues: ValidationIssue[] = []
    const atomMap = new Map(atoms.map((atom) => [atom.id, atom]))
    bonds.forEach((bond) => {
      const atomA = atomMap.get(bond.a1)
      const atomB = atomMap.get(bond.a2)
      if (!atomA || !atomB) return
      const distance = this.distance(atomA.position, atomB.position)
      const expected = approximateBondLength(atomA.element, atomB.element, bond.order)
      if (Math.abs(distance - expected) > 0.25) {
        issues.push({
          id: generateId(),
          type: 'bond_length',
          severity: 'warning',
          message: `Unusual ${atomA.element}-${atomB.element} bond length`,
          atoms: [atomA.id, atomB.id],
          bonds: [bond.id],
          expectedValue: parseFloat(expected.toFixed(2)),
          actualValue: parseFloat(distance.toFixed(2)),
          recommendation: 'Consider optimizing bond geometry or adjusting coordinates.',
        })
      }
    })
    return issues
  }

  private checkAngleStrain(
    atoms: SimpleAtom[],
    bonds: SimpleBond[]
  ): ValidationIssue[] {
    const issues: ValidationIssue[] = []
    const adjacency = new Map<string, SimpleBond[]>()
    bonds.forEach((bond) => {
      adjacency.set(bond.a1, [...(adjacency.get(bond.a1) ?? []), bond])
      adjacency.set(bond.a2, [...(adjacency.get(bond.a2) ?? []), bond])
    })
    const atomMap = new Map(atoms.map((atom) => [atom.id, atom]))
    atoms.forEach((atom) => {
      const connected = adjacency.get(atom.id) ?? []
      if (connected.length < 2) return
      for (let i = 0; i < connected.length; i++) {
        for (let j = i + 1; j < connected.length; j++) {
          const neighborA = this.otherAtomId(connected[i], atom.id)
          const neighborB = this.otherAtomId(connected[j], atom.id)
          const atomA = atomMap.get(neighborA)
          const atomB = atomMap.get(neighborB)
          if (!atomA || !atomB) continue
          const angle = this.bondAngle(atomA.position, atom.position, atomB.position)
          if (angle < 60 || angle > 150) {
            issues.push({
              id: generateId(),
              type: 'angle_strain',
              severity: 'warning',
              message: `Angle strain detected around ${atom.element}`,
              atoms: [atom.id, neighborA, neighborB],
              expectedValue: 109.5,
              actualValue: parseFloat(angle.toFixed(1)),
              recommendation: 'Adjust substituent positions or run geometry optimization.',
            })
          }
        }
      }
    })
    return issues
  }

  private generateSuggestions(issues: ValidationIssue[]): string[] {
    if (!issues.length) return ['Structure looks stable.']
    const suggestions = new Set<string>()
    issues.forEach((issue) => {
      if (issue.type === 'steric_clash') {
        suggestions.add('Resolve steric clashes by increasing distances or optimizing geometry.')
      }
      if (issue.type === 'bond_length') {
        suggestions.add('Adjust unusual bond lengths to typical values.')
      }
      if (issue.type === 'angle_strain') {
        suggestions.add('Relieve angle strain by repositioning substituents.')
      }
    })
    return Array.from(suggestions)
  }

  private distance(a: [number, number, number], b: [number, number, number]) {
    return Math.sqrt(
      (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2
    )
  }

  private otherAtomId(bond: SimpleBond, atomId: string) {
    return bond.a1 === atomId ? bond.a2 : bond.a1
  }

  private bondAngle(
    a: [number, number, number],
    b: [number, number, number],
    c: [number, number, number]
  ) {
    const ab = [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
    const cb = [c[0] - b[0], c[1] - b[1], c[2] - b[2]]
    const dot = ab[0] * cb[0] + ab[1] * cb[1] + ab[2] * cb[2]
    const magAB = Math.sqrt(ab[0] ** 2 + ab[1] ** 2 + ab[2] ** 2)
    const magCB = Math.sqrt(cb[0] ** 2 + cb[1] ** 2 + cb[2] ** 2)
    const cos = dot / (magAB * magCB || 1)
    const angle = (Math.acos(Math.min(1, Math.max(-1, cos))) * 180) / Math.PI
    return angle
  }
}

export class StructureOptimizer {
  async optimize(molecule: MoleculeGraph): Promise<MoleculeGraph> {
    const clone = molecule.clone()
    ForceField.optimizeGeometry(clone, 25, 0.01)
    return clone
  }
}


