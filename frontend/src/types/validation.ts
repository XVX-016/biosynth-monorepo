export type ValidationIssueType =
  | 'steric_clash'
  | 'bond_length'
  | 'angle_strain'
  | 'chirality'
  | 'charge'

export interface ValidationIssue {
  id: string
  type: ValidationIssueType
  severity: 'warning' | 'error'
  message: string
  atoms: string[]
  bonds?: string[]
  expectedValue?: number
  actualValue?: number
  recommendation?: string
}

export interface AutoFix {
  id: string
  label: string
  description: string
  appliesTo: string[]
}

export interface ValidationResult {
  issues: ValidationIssue[]
  score: number
  suggestions: string[]
  autoFixes: AutoFix[]
  timestamp: string
}

export interface ValidationState {
  status: 'idle' | 'running' | 'complete' | 'error'
  result: ValidationResult | null
  error?: string | null
  lastValidated?: string
}



