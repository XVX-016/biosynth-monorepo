/**
 * Molecule validation engine
 * 
 * Phase 5: Complete validation implementation
 * 
 * Validates molecular structure for:
 * - Valence violations
 * - Bond order issues
 * - Structural problems
 * - Ring strain
 * - Overlaps
 * - Disconnected fragments
 */

import type { ValidationResult, ValidationError, ValidationErrorCode } from '../types'
import type { Molecule } from '../Molecule'
import { ELEMENT_DATA } from '../constants'

export function validateMolecule(molecule: Molecule): ValidationResult {
  const errors: ValidationError[] = []
  const warnings: ValidationError[] = []

  if (molecule.isEmpty()) {
    return {
      valid: true,
      errors: [],
      warnings: [],
    }
  }

  // 1. Valence checking
  const valenceErrors = checkValence(molecule)
  errors.push(...valenceErrors)

  // 2. Bond order validation
  const bondErrors = checkBondOrders(molecule)
  errors.push(...bondErrors)

  // 3. Overlap detection
  const overlapWarnings = checkOverlaps(molecule)
  warnings.push(...overlapWarnings)

  // 4. Disconnected fragment detection
  const fragmentWarnings = checkFragments(molecule)
  warnings.push(...fragmentWarnings)

  // 5. Ring strain detection (basic)
  const strainWarnings = checkRingStrain(molecule)
  warnings.push(...strainWarnings)

  return {
    valid: errors.length === 0,
    errors,
    warnings,
  }
}

/**
 * Check valence violations
 */
function checkValence(molecule: Molecule): ValidationError[] {
  const errors: ValidationError[] = []
  const atoms = molecule.getAtoms()

  atoms.forEach(atom => {
    const bonds = molecule.getBondsForAtom(atom.id)
    const totalBondOrder = bonds.reduce((sum, bond) => sum + bond.order, 0)
    const elementInfo = ELEMENT_DATA[atom.element]

    if (!elementInfo) {
      errors.push({
        code: 'INVALID_ELEMENT',
        message: `Unknown element: ${atom.element}`,
        atomId: atom.id,
        details: { element: atom.element },
      })
      return
    }

    const maxValence = elementInfo.maxValence
    const currentValence = totalBondOrder

    if (currentValence > maxValence) {
      errors.push({
        code: 'VALENCE_EXCEEDED',
        message: `${atom.element} atom exceeds maximum valence (${currentValence} > ${maxValence})`,
        atomId: atom.id,
        details: {
          element: atom.element,
          currentValence,
          maxValence,
          allowed: maxValence,
          got: currentValence,
        },
      })
    }

    // Check for underfulfilled valence (warning, not error)
    const commonValences = elementInfo.commonValence
    if (commonValences.length > 0) {
      const minCommonValence = Math.min(...commonValences)
      if (currentValence < minCommonValence && currentValence > 0) {
        // Only warn if not zero (zero is valid for noble gases)
        errors.push({
          code: 'VALENCE_UNDERFULFILLED',
          message: `${atom.element} atom may have underfulfilled valence (${currentValence} < ${minCommonValence})`,
          atomId: atom.id,
          details: {
            element: atom.element,
            currentValence,
            minCommonValence,
          },
        })
      }
    }
  })

  return errors
}

/**
 * Check bond order validity
 */
function checkBondOrders(molecule: Molecule): ValidationError[] {
  const errors: ValidationError[] = []
  const bonds = molecule.getBonds()

  bonds.forEach(bond => {
    const atom1 = molecule.getAtom(bond.atom1)
    const atom2 = molecule.getAtom(bond.atom2)

    if (!atom1 || !atom2) {
      errors.push({
        code: 'INVALID_BOND_ORDER',
        message: 'Bond connects non-existent atoms',
        bondId: bond.id,
        details: { atom1: bond.atom1, atom2: bond.atom2 },
      })
      return
    }

    // Check bond order is valid (1, 2, 3, or 1.5 for aromatic)
    if (![1, 2, 3, 1.5].includes(bond.order)) {
      errors.push({
        code: 'INVALID_BOND_ORDER',
        message: `Invalid bond order: ${bond.order}`,
        bondId: bond.id,
        details: { order: bond.order },
      })
    }

    // Check if bond order is too high for the elements
    const element1Info = ELEMENT_DATA[atom1.element]
    const element2Info = ELEMENT_DATA[atom2.element]

    if (element1Info && element2Info) {
      // Triple bonds are rare, mostly for C, N
      if (bond.order === 3) {
        const canHaveTriple = ['C', 'N'].includes(atom1.element) && ['C', 'N'].includes(atom2.element)
        if (!canHaveTriple) {
          errors.push({
            code: 'INVALID_BOND_ORDER',
            message: `${atom1.element}-${atom2.element} triple bond is unusual`,
            bondId: bond.id,
            details: { order: bond.order, element1: atom1.element, element2: atom2.element },
          })
        }
      }
    }
  })

  return errors
}

/**
 * Check for atom overlaps (atoms too close together)
 */
function checkOverlaps(molecule: Molecule): ValidationError[] {
  const warnings: ValidationError[] = []
  const atoms = molecule.getAtoms()
  const OVERLAP_THRESHOLD = 0.5 // Angstroms

  for (let i = 0; i < atoms.length; i++) {
    for (let j = i + 1; j < atoms.length; j++) {
      const atom1 = atoms[i]
      const atom2 = atoms[j]

      // Check if atoms are bonded (bonded atoms can be close)
      const bonds = molecule.getBonds()
      const areBonded = bonds.some(bond =>
        (bond.atom1 === atom1.id && bond.atom2 === atom2.id) ||
        (bond.atom1 === atom2.id && bond.atom2 === atom1.id)
      )

      if (areBonded) continue

      // Calculate distance
      const [x1, y1, z1] = atom1.position
      const [x2, y2, z2] = atom2.position
      const distance = Math.sqrt(
        Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2) + Math.pow(z2 - z1, 2)
      )

      if (distance < OVERLAP_THRESHOLD) {
        warnings.push({
          code: 'BOND_OVERLAP',
          message: `Atoms ${atom1.element} and ${atom2.element} are very close (${distance.toFixed(2)} Å)`,
          atomId: atom1.id,
          details: {
            atom1Id: atom1.id,
            atom2Id: atom2.id,
            distance,
            threshold: OVERLAP_THRESHOLD,
          },
        })
      }
    }
  }

  return warnings
}

/**
 * Check for disconnected fragments
 */
function checkFragments(molecule: Molecule): ValidationError[] {
  const warnings: ValidationError[] = []

  if (molecule.isEmpty()) {
    return warnings
  }

  const atoms = molecule.getAtoms()
  const visited = new Set<string>()
  const fragments: string[][] = []

  // Find all connected components
  atoms.forEach(atom => {
    if (visited.has(atom.id)) return

    const fragment: string[] = []
    const queue: string[] = [atom.id]

    while (queue.length > 0) {
      const currentId = queue.shift()!
      if (visited.has(currentId)) continue

      visited.add(currentId)
      fragment.push(currentId)

      // Add neighbors
      const neighbors = molecule.getNeighbors(currentId)
      neighbors.forEach(neighbor => {
        if (!visited.has(neighbor.id)) {
          queue.push(neighbor.id)
        }
      })
    }

    fragments.push(fragment)
  })

  // Warn if multiple fragments
  if (fragments.length > 1) {
    warnings.push({
      code: 'DISCONNECTED_FRAGMENT',
      message: `Molecule has ${fragments.length} disconnected fragments`,
      details: {
        fragmentCount: fragments.length,
        fragmentSizes: fragments.map(f => f.length),
      },
    })
  }

  return warnings
}

/**
 * Check for ring strain (basic detection)
 */
function checkRingStrain(molecule: Molecule): ValidationError[] {
  const warnings: ValidationError[] = []

  // Basic ring detection: find cycles
  const atoms = molecule.getAtoms()
  const cycles = findCycles(molecule)

  cycles.forEach(cycle => {
    if (cycle.length < 3) {
      warnings.push({
        code: 'RING_STRAIN',
        message: `Ring with ${cycle.length} atoms is impossible`,
        details: { cycleSize: cycle.length, atoms: cycle },
      })
    } else if (cycle.length === 3) {
      warnings.push({
        code: 'RING_STRAIN',
        message: '3-membered ring has high strain',
        details: { cycleSize: 3, atoms: cycle },
      })
    } else if (cycle.length === 4) {
      warnings.push({
        code: 'RING_STRAIN',
        message: '4-membered ring has moderate strain',
        details: { cycleSize: 4, atoms: cycle },
      })
    }
  })

  return warnings
}

/**
 * Find cycles in molecule (simple DFS)
 */
function findCycles(molecule: Molecule): string[][] {
  const cycles: string[][] = []
  const atoms = molecule.getAtoms()
  const visited = new Set<string>()

  atoms.forEach(startAtom => {
    if (visited.has(startAtom.id)) return

    const path: string[] = []
    const findCycle = (currentId: string, parentId: string | null) => {
      if (visited.has(currentId)) {
        // Found a cycle
        const cycleStart = path.indexOf(currentId)
        if (cycleStart !== -1) {
          const cycle = path.slice(cycleStart)
          if (cycle.length >= 3) {
            cycles.push([...cycle, currentId])
          }
        }
        return
      }

      visited.add(currentId)
      path.push(currentId)

      const neighbors = molecule.getNeighbors(currentId)
      neighbors.forEach(neighbor => {
        if (neighbor.id !== parentId) {
          findCycle(neighbor.id, currentId)
        }
      })

      path.pop()
    }

    findCycle(startAtom.id, null)
  })

  // Remove duplicate cycles
  const uniqueCycles: string[][] = []
  cycles.forEach(cycle => {
    const sorted = [...cycle].sort().join(',')
    if (!uniqueCycles.some(c => [...c].sort().join(',') === sorted)) {
      uniqueCycles.push(cycle)
    }
  })

  return uniqueCycles
}

/**
 * Validate atom can accept bond
 */
export function canAtomAcceptBond(molecule: Molecule, atomId: string, bondOrder: number = 1): boolean {
  const atom = molecule.getAtom(atomId)
  if (!atom) return false

  const bonds = molecule.getBondsForAtom(atomId)
  const currentValence = bonds.reduce((sum, bond) => sum + bond.order, 0)
  const elementInfo = ELEMENT_DATA[atom.element]

  if (!elementInfo) return false

  return (currentValence + bondOrder) <= elementInfo.maxValence
}

/**
 * Validate bond can be created
 */
export function canCreateBond(
  molecule: Molecule,
  atom1Id: string,
  atom2Id: string,
  order: number = 1
): { valid: boolean; error?: ValidationError } {
  // Check atoms exist
  const atom1 = molecule.getAtom(atom1Id)
  const atom2 = molecule.getAtom(atom2Id)

  if (!atom1 || !atom2) {
    return {
      valid: false,
      error: {
        code: 'INVALID_BOND_ORDER',
        message: 'One or both atoms do not exist',
        details: { atom1Id, atom2Id },
      },
    }
  }

  // Check self-bond
  if (atom1Id === atom2Id) {
    return {
      valid: false,
      error: {
        code: 'INVALID_BOND_ORDER',
        message: 'Cannot bond atom to itself',
        details: { atom1Id, atom2Id },
      },
    }
  }

  // Check if bond already exists
  const existingBonds = molecule.getBonds()
  for (const bond of existingBonds) {
    if (
      (bond.atom1 === atom1Id && bond.atom2 === atom2Id) ||
      (bond.atom1 === atom2Id && bond.atom2 === atom1Id)
    ) {
      return {
        valid: false,
        error: {
          code: 'INVALID_BOND_ORDER',
          message: 'Bond already exists',
          bondId: bond.id,
          details: { atom1Id, atom2Id },
        },
      }
    }
  }

  // Check valence
  if (!canAtomAcceptBond(molecule, atom1Id, order)) {
    return {
      valid: false,
      error: {
        code: 'VALENCE_EXCEEDED',
        message: `${atom1.element} atom cannot accept bond of order ${order}`,
        atomId: atom1Id,
        details: { element: atom1.element, bondOrder: order },
      },
    }
  }

  if (!canAtomAcceptBond(molecule, atom2Id, order)) {
    return {
      valid: false,
      error: {
        code: 'VALENCE_EXCEEDED',
        message: `${atom2.element} atom cannot accept bond of order ${order}`,
        atomId: atom2Id,
        details: { element: atom2.element, bondOrder: order },
      },
    }
  }

  return { valid: true }
}
