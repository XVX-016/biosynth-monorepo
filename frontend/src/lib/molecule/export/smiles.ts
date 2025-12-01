/**
 * SMILES export functions
 * 
 * Phase 6: Integration with RDKit Backend
 * 
 * Converts molecule to SMILES string via backend API.
 */

import type { Molecule } from '../Molecule'

/**
 * Convert molecule to SMILES string
 * 
 * Uses backend API for accurate SMILES generation with RDKit.
 */
export async function toSMILES(molecule: Molecule, canonicalize: boolean = true): Promise<string> {
  if (molecule.isEmpty()) {
    return ''
  }

  // Convert molecule to JSON format for backend
  const moleculeData = moleculeToBackendFormat(molecule)

  try {
    const response = await fetch('/api/molecule/to-smiles', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        molecule: moleculeData,
        canonicalize,
      }),
    })

    if (!response.ok) {
      throw new Error(`Failed to generate SMILES: ${response.statusText}`)
    }

    const result = await response.json()
    return result.smiles || ''
  } catch (error) {
    console.error('Error generating SMILES:', error)
    // Fallback: generate basic SMILES from structure
    return generateBasicSMILES(molecule)
  }
}

/**
 * Convert molecule to MolBlock format
 */
export async function toMolBlock(molecule: Molecule): Promise<string> {
  if (molecule.isEmpty()) {
    return ''
  }

  const moleculeData = moleculeToBackendFormat(molecule)

  try {
    const response = await fetch('/api/molecule/to-molblock', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        molecule: moleculeData,
      }),
    })

    if (!response.ok) {
      throw new Error(`Failed to generate MolBlock: ${response.statusText}`)
    }

    const result = await response.json()
    return result.molblock || ''
  } catch (error) {
    console.error('Error generating MolBlock:', error)
    throw error
  }
}

/**
 * Normalize hydrogens in molecule
 */
export async function normalizeHydrogens(molecule: Molecule): Promise<Molecule> {
  if (molecule.isEmpty()) {
    return molecule
  }

  const moleculeData = moleculeToBackendFormat(molecule)

  try {
    const response = await fetch('/api/molecule/normalize-hydrogens', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        molecule: moleculeData,
      }),
    })

    if (!response.ok) {
      throw new Error(`Failed to normalize hydrogens: ${response.statusText}`)
    }

    const result = await response.json()
    // Convert back to Molecule
    return moleculeFromBackendFormat(result.molecule)
  } catch (error) {
    console.error('Error normalizing hydrogens:', error)
    return molecule // Return original on error
  }
}

/**
 * Validate molecule with RDKit backend
 */
export async function validateWithRDKit(molecule: Molecule): Promise<{
  valid: boolean
  sanitized_smiles?: string
  molblock?: string
  errors: string[]
}> {
  if (molecule.isEmpty()) {
    return { valid: true, errors: [] }
  }

  const moleculeData = moleculeToBackendFormat(molecule)

  try {
    const response = await fetch('/api/molecule/validate', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        molecule: moleculeData,
      }),
    })

    if (!response.ok) {
      throw new Error(`Validation failed: ${response.statusText}`)
    }

    return await response.json()
  } catch (error) {
    console.error('Error validating molecule:', error)
    return {
      valid: false,
      errors: [error instanceof Error ? error.message : 'Validation failed'],
    }
  }
}

/**
 * Convert molecule to backend format
 */
function moleculeToBackendFormat(molecule: Molecule): any {
  const atoms = molecule.getAtoms().map(atom => ({
    id: atom.id,
    element: atom.element,
    position: atom.position,
    charge: atom.charge,
    formalCharge: atom.formalCharge,
  }))

  const bonds = molecule.getBonds().map(bond => ({
    id: bond.id,
    atom1: bond.atom1,
    atom2: bond.atom2,
    order: bond.order,
    type: bond.type,
  }))

  return {
    atoms,
    bonds,
    metadata: molecule.getMetadata(),
  }
}

/**
 * Convert backend format to Molecule
 */
function moleculeFromBackendFormat(data: any): Molecule {
  const { Molecule } = require('../Molecule')
  const molecule = new Molecule()

  // Add atoms
  data.atoms?.forEach((atom: any) => {
    molecule.addAtom({
      id: atom.id,
      element: atom.element,
      position: atom.position,
      charge: atom.charge,
      formalCharge: atom.formalCharge,
    })
  })

  // Add bonds
  data.bonds?.forEach((bond: any) => {
    molecule.addBond({
      id: bond.id,
      atom1: bond.atom1,
      atom2: bond.atom2,
      order: bond.order,
      type: bond.type,
    })
  })

  if (data.metadata) {
    molecule.setMetadata(data.metadata)
  }

  return molecule
}

/**
 * Generate basic SMILES (fallback, not accurate)
 * 
 * This is a simple fallback when backend is unavailable.
 * For accurate SMILES, use toSMILES() with backend.
 */
function generateBasicSMILES(molecule: Molecule): string {
  // Very basic implementation - just concatenate elements
  // This should never be used in production, only as fallback
  const atoms = molecule.getAtoms()
  if (atoms.length === 0) return ''

  const elements = atoms.map(atom => atom.element).join('')
  return `[${elements}]` // Placeholder format
}

