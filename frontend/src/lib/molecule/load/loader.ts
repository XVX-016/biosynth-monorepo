/**
 * Molecule loader functions
 * 
 * Phase 13: Save / Load / Export System
 * 
 * Handles loading molecules from various formats.
 */

import type { Molecule } from '../Molecule'
import { fromJSON } from '../export/json'

// Import moleculeFromBackendFormat from smiles.ts
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
 * Load molecule from SMILES string
 */
export async function loadFromSMILES(smiles: string): Promise<Molecule> {
  if (!smiles || !smiles.trim()) {
    throw new Error('Empty SMILES string')
  }

  try {
    // Use backend API to convert SMILES to molecule
    const response = await fetch('/api/molecule/from-smiles', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ smiles: smiles.trim() }),
    })

    if (!response.ok) {
      throw new Error(`Failed to load SMILES: ${response.statusText}`)
    }

    const result = await response.json()
    return moleculeFromBackendFormat(result.molecule)
  } catch (error) {
    console.error('Error loading from SMILES:', error)
    throw error
  }
}

/**
 * Load molecule from MolBlock string
 */
export async function loadFromMolBlock(molblock: string): Promise<Molecule> {
  if (!molblock || !molblock.trim()) {
    throw new Error('Empty MolBlock string')
  }

  try {
    // Use backend API to convert MolBlock to molecule
    const response = await fetch('/api/molecule/from-molblock', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ molblock: molblock.trim() }),
    })

    if (!response.ok) {
      throw new Error(`Failed to load MolBlock: ${response.statusText}`)
    }

    const result = await response.json()
    return moleculeFromBackendFormat(result.molecule)
  } catch (error) {
    console.error('Error loading from MolBlock:', error)
    throw error
  }
}

/**
 * Load molecule from JSON string
 */
export function loadFromJSON(json: string): Molecule {
  return fromJSON(json)
}

/**
 * Load molecule from file
 */
export async function loadFromFile(file: File): Promise<Molecule> {
  const text = await file.text()
  const extension = file.name.split('.').pop()?.toLowerCase()

  switch (extension) {
    case 'json':
      return loadFromJSON(text)
    case 'mol':
    case 'sdf':
      return loadFromMolBlock(text)
    case 'smi':
    case 'smiles':
      return loadFromSMILES(text)
    default:
      // Try to detect format
      if (text.trim().startsWith('{')) {
        // JSON
        return loadFromJSON(text)
      } else if (text.includes('V2000') || text.includes('V3000')) {
        // MolBlock
        return loadFromMolBlock(text)
      } else {
        // Assume SMILES
        return loadFromSMILES(text)
      }
  }
}

/**
 * Load molecule from drag-and-drop data
 */
export async function loadFromDropEvent(event: DragEvent): Promise<Molecule | null> {
  event.preventDefault()
  
  const files = event.dataTransfer?.files
  if (!files || files.length === 0) {
    return null
  }

  const file = files[0]
  return loadFromFile(file)
}

