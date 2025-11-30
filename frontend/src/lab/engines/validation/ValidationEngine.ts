/**
 * ValidationEngine - Main validation engine
 */

import { ValidationResult, ValidationError } from './Validation.types';
import { validateValence } from './BasicValenceValidator';
import { runAdvancedValidations } from './AdvancedValidator';
import { sanitize } from './Sanitizer';
import type { MoleculeState } from '../MoleculeStateEngine';

export class ValidationEngine {
  /**
   * Validate molecule state
   */
  validate(state: MoleculeState): ValidationResult {
    const errors: ValidationError[] = [];
    const warnings: string[] = [];

    // Basic valence validation
    const valenceErrors = validateValence(state);
    errors.push(...valenceErrors);

    // Advanced validations
    const advancedErrors = runAdvancedValidations(state);
    errors.push(...advancedErrors);

    // Check for disconnected atoms
    const atomIds = new Set(Array.from(state.atoms.keys()));
    const connectedAtoms = new Set<string>();
    
    state.bonds.forEach((bond) => {
      connectedAtoms.add(bond.atoms[0]);
      connectedAtoms.add(bond.atoms[1]);
    });

    atomIds.forEach((atomId) => {
      if (!connectedAtoms.has(atomId) && state.atoms.size > 1) {
        warnings.push(`Atom ${atomId} is disconnected from the molecule`);
      }
    });

    // Check for duplicate bonds
    const bondKeys = new Set<string>();
    state.bonds.forEach((bond) => {
      const key = [bond.atoms[0], bond.atoms[1]].sort().join('-');
      if (bondKeys.has(key)) {
        errors.push({
          type: 'duplicate_bond',
          message: 'Duplicate bond detected',
          bondId: bond.id,
        });
      }
      bondKeys.add(key);
    });

    return {
      errors,
      warnings,
      valid: errors.length === 0,
    };
  }

  /**
   * Sanitize molecule state
   */
  sanitize(state: MoleculeState): MoleculeState {
    return sanitize(state);
  }

  /**
   * Validate and sanitize
   */
  validateAndSanitize(state: MoleculeState): ValidationResult {
    const validation = this.validate(state);
    if (validation.valid) {
      validation.sanitizedState = this.sanitize(state);
    }
    return validation;
  }
}

