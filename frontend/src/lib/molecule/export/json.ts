/**
 * JSON export/import functions
 * 
 * Phase 13: Save / Load / Export System
 * 
 * Handles molecule serialization to/from JSON format.
 */

import type { Molecule } from '../Molecule'
import type { MoleculeState } from '../types'

/**
 * Export molecule to JSON
 */
export function toJSON(molecule: Molecule): string {
  const state: MoleculeState = molecule.toState()
  return JSON.stringify(state, null, 2)
}

/**
 * Import molecule from JSON
 */
export function fromJSON(json: string): Molecule {
  const state: MoleculeState = JSON.parse(json)
  const { Molecule } = require('../Molecule')
  return new Molecule(state)
}

/**
 * Export molecule to JSON file
 */
export function exportMoleculeAsJSON(molecule: Molecule, filename: string = 'molecule.json'): void {
  const json = toJSON(molecule)
  const blob = new Blob([json], { type: 'application/json' })
  const url = URL.createObjectURL(blob)
  const link = document.createElement('a')
  link.href = url
  link.download = filename
  document.body.appendChild(link)
  link.click()
  document.body.removeChild(link)
  URL.revokeObjectURL(url)
}

