/**
 * Validation Engine - Types
 */

export interface ValidationError {
  type: string;
  message: string;
  atomId?: string;
  bondId?: string;
}

export interface ValidationResult {
  errors: ValidationError[];
  warnings: string[];
  sanitizedState?: any;
  valid: boolean;
}

