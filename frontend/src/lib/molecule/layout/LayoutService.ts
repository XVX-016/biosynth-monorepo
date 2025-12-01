/**
 * Layout Service
 * 
 * Phase 8: 2D Layout Generation
 * 
 * Provides 2D coordinate generation using RDKit backend:
 * - Auto-layout on import
 * - Straighten chains
 * - Symmetrize rings
 * - Spacing correction
 */

import type { Molecule } from '../Molecule'

// Import conversion functions (avoid circular dependency)
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

function moleculeFromBackendFormat(data: any): Molecule {
  const { Molecule } = require('../Molecule')
  const molecule = new Molecule()

  // Add atoms
  data.atoms?.forEach((atom: any) => {
    molecule.addAtom({
      id: atom.id,
      element: atom.element,
      position: atom.position || [0, 0, 0],
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
      order: bond.order || 1,
      type: bond.type,
    })
  })

  if (data.metadata) {
    molecule.setMetadata(data.metadata)
  }

  return molecule
}

export interface LayoutOptions {
  method?: 'coordgen' | 'rdkit'
  spacing?: number
}

/**
 * Generate 2D coordinates for molecule
 */
export async function generate2DLayout(
  molecule: Molecule,
  options: LayoutOptions = {}
): Promise<Molecule> {
  if (molecule.isEmpty()) {
    return molecule
  }

  const moleculeData = moleculeToBackendFormat(molecule)

  try {
    const response = await fetch('/api/molecule/generate-2d-layout', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        molecule: moleculeData,
        method: options.method || 'coordgen',
        spacing: options.spacing || 1.5,
      }),
    })

    if (!response.ok) {
      throw new Error(`Failed to generate 2D layout: ${response.statusText}`)
    }

    const result = await response.json()
    return moleculeFromBackendFormat(result.molecule)
  } catch (error) {
    console.error('Error generating 2D layout:', error)
    // Fallback: return original molecule
    return molecule
  }
}

/**
 * Generate 2D layout from SMILES
 */
export async function generate2DLayoutFromSMILES(
  smiles: string,
  options: LayoutOptions = {}
): Promise<Molecule> {
  try {
    const response = await fetch('/api/molecule/generate-2d-layout', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        smiles,
        method: options.method || 'coordgen',
        spacing: options.spacing || 1.5,
      }),
    })

    if (!response.ok) {
      throw new Error(`Failed to generate 2D layout: ${response.statusText}`)
    }

    const result = await response.json()
    return moleculeFromBackendFormat(result.molecule)
  } catch (error) {
    console.error('Error generating 2D layout from SMILES:', error)
    throw error
  }
}

