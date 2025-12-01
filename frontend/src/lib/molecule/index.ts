/**
 * Molecule Editor Library
 * 
 * Centralized molecule editing, validation, and operations.
 * 
 * This is the single source of truth for molecule manipulation.
 */

// Core types
export type {
  Atom,
  Bond,
  MoleculeState,
  EditorTool,
  ValidationError,
  ValidationErrorCode,
  ValidationResult,
  ElementInfo,
  ExportOptions,
} from './types'

// Core classes
export { AtomImpl } from './Atom'
export { BondImpl } from './Bond'
export { Molecule } from './Molecule'

// Constants
export { ELEMENT_DATA, COMMON_ELEMENTS, DEFAULT_BOND_ORDER, VALID_BOND_ORDERS } from './constants'

// Validation
export { validateMolecule, canAtomAcceptBond, canCreateBond } from './validation/Validator'
export type { ValidationResult, ValidationError, ValidationErrorCode } from './types'

// Export functions (RDKit backend integration)
export { toSMILES, toMolBlock, normalizeHydrogens, validateWithRDKit } from './export'

// Prediction service (ML pipeline)
export { predictionService, PredictionService } from './prediction'
export type { PredictionResult, AttentionData, PredictionRequest } from './prediction'

// Layout service (2D coordinate generation)
export { generate2DLayout, generate2DLayoutFromSMILES } from './layout'
export type { LayoutOptions } from './layout'

