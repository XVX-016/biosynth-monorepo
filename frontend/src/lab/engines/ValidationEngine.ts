/**
 * ValidationEngine - Chemical structure validation
 * 
 * Provides both local quick validation and deep validation via backend RDKit
 */

import { moleculeEngine, type Atom, type Bond } from './MoleculeStateEngine'

export interface ValidationIssue {
  type: 'error' | 'warning' | 'info'
  message: string
  atomId?: string
  bondId?: string
}

export interface ValidationResult {
  valid: boolean
  issues: ValidationIssue[]
}

/**
 * Valence rules for common elements
 */
const VALENCE_RULES: Record<string, number> = {
  H: 1,
  C: 4,
  N: 3,
  O: 2,
  F: 1,
  Cl: 1,
  Br: 1,
  I: 1,
  S: 6,
  P: 5,
  B: 3,
  Si: 4,
}

export class ValidationEngine {
  /**
   * Quick local validation (frontend-only checks)
   */
  static localCheck(): ValidationResult {
    const issues: ValidationIssue[] = []

    // Check for disconnected atoms
    const atoms = moleculeEngine.getAllAtoms()
    const bonds = moleculeEngine.getAllBonds()
    
    atoms.forEach(atom => {
      const connectedBonds = moleculeEngine.getBondsForAtom(atom.id)
      
      // Check valence
      const expectedValence = VALENCE_RULES[atom.element]
      if (expectedValence !== undefined) {
        const actualValence = connectedBonds.reduce((sum, bond) => sum + bond.order, 0)
        if (actualValence > expectedValence) {
          issues.push({
            type: 'error',
            message: `${atom.element} atom has valence ${actualValence}, expected max ${expectedValence}`,
            atomId: atom.id
          })
        } else if (actualValence < expectedValence && atom.element !== 'H') {
          issues.push({
            type: 'warning',
            message: `${atom.element} atom has incomplete valence (${actualValence}/${expectedValence})`,
            atomId: atom.id
          })
        }
      }

      // Check for disconnected atoms (no bonds)
      if (connectedBonds.length === 0 && atoms.length > 1) {
        issues.push({
          type: 'warning',
          message: `${atom.element} atom is disconnected`,
          atomId: atom.id
        })
      }
    })

    // Check for duplicate bonds
    const bondPairs = new Set<string>()
    bonds.forEach(bond => {
      const pair = [bond.atoms[0], bond.atoms[1]].sort().join('-')
      if (bondPairs.has(pair)) {
        issues.push({
          type: 'error',
          message: 'Duplicate bond detected',
          bondId: bond.id
        })
      }
      bondPairs.add(pair)
    })

    // Check for self-bonds
    bonds.forEach(bond => {
      if (bond.atoms[0] === bond.atoms[1]) {
        issues.push({
          type: 'error',
          message: 'Atom cannot bond to itself',
          bondId: bond.id
        })
      }
    })

    return {
      valid: issues.filter(i => i.type === 'error').length === 0,
      issues
    }
  }

  /**
   * Deep validation using backend RDKit
   */
  static async deepValidate(smiles: string): Promise<ValidationResult> {
    try {
      const response = await fetch('/api/validate', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ smiles }),
      })

      if (!response.ok) {
        throw new Error(`Validation failed: ${response.statusText}`)
      }

      const data = await response.json()
      
      return {
        valid: data.valid || false,
        issues: (data.issues || []).map((issue: any) => ({
          type: issue.type || 'warning',
          message: issue.message || 'Unknown issue',
          atomId: issue.atomId,
          bondId: issue.bondId,
        }))
      }
    } catch (error) {
      console.error('Validation error:', error)
      return {
        valid: false,
        issues: [{
          type: 'error',
          message: error instanceof Error ? error.message : 'Validation service unavailable'
        }]
      }
    }
  }

  /**
   * Validate current molecule state
   */
  static async validateCurrent(): Promise<ValidationResult> {
    // First do local check
    const localResult = this.localCheck()
    
    // If local validation passes, try deep validation
    if (localResult.valid) {
      try {
        const smiles = await moleculeEngine.serializeSMILES()
        if (smiles && smiles !== 'C') {
          const deepResult = await this.deepValidate(smiles)
          return {
            valid: deepResult.valid,
            issues: [...localResult.issues, ...deepResult.issues]
          }
        }
      } catch (error) {
        console.warn('Deep validation failed, using local result:', error)
      }
    }
    
    return localResult
  }
}

