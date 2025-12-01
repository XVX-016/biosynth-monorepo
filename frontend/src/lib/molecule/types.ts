/**
 * Core molecule data types
 * 
 * Centralized definitions for atoms, bonds, and molecule graph structure.
 * This is the single source of truth for molecular data structures.
 */

export interface Atom {
  id: string
  element: string
  position: [number, number, number]  // x, y, z coordinates
  charge?: number
  formalCharge?: number
  hybridization?: string
  aromatic?: boolean
  inRing?: boolean
  valence?: number
  maxValence?: number
  metadata?: Record<string, any>
}

export interface Bond {
  id: string
  atom1: string  // atom id
  atom2: string  // atom id
  order: number  // 1, 2, 3, or aromatic
  type?: 'single' | 'double' | 'triple' | 'aromatic'
  stereo?: string
  metadata?: Record<string, any>
}

export interface MoleculeState {
  atoms: Map<string, Atom>
  bonds: Map<string, Bond>
  metadata?: {
    name?: string
    formula?: string
    smiles?: string
    molblock?: string
    [key: string]: any
  }
}

/**
 * Tool modes for the molecule editor
 */
export type EditorTool = 
  | 'select'      // Select atoms/bonds
  | 'add-atom'    // Place new atoms
  | 'bond'        // Create bonds between atoms
  | 'delete'      // Delete atoms/bonds
  | 'move'        // Move atoms
  | 'inspect'     // Inspect properties

/**
 * Validation error types
 */
export type ValidationErrorCode =
  | 'VALENCE_EXCEEDED'
  | 'VALENCE_UNDERFULFILLED'
  | 'INVALID_BOND_ORDER'
  | 'BOND_OVERLAP'
  | 'DISCONNECTED_FRAGMENT'
  | 'RING_STRAIN'
  | 'INVALID_ELEMENT'
  | 'CHARGE_IMBALANCE'

export interface ValidationError {
  code: ValidationErrorCode
  message: string
  atomId?: string
  bondId?: string
  details?: Record<string, any>
}

export interface ValidationResult {
  valid: boolean
  errors: ValidationError[]
  warnings: ValidationError[]
}

/**
 * Element information
 */
export interface ElementInfo {
  symbol: string
  atomicNumber: number
  name: string
  maxValence: number
  commonValence: number[]
  electronegativity?: number
  radius?: number
}

/**
 * Export formats
 */
export interface ExportOptions {
  format: 'smiles' | 'molblock' | 'json' | 'coordinates'
  includeHydrogens?: boolean
  canonicalize?: boolean
}

