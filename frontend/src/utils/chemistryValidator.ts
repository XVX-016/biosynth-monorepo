/**
 * Chemistry Validator
 * Validates that prompts are chemistry-related
 */

const CHEMISTRY_KEYWORDS = [
  "molecule",
  "reaction",
  "bond",
  "atom",
  "smiles",
  "synthesis",
  "chemistry",
  "pKa",
  "solubility",
  "mechanism",
  "functional group",
  "organic",
  "inorganic",
  "compound",
  "chemical",
  "molecular",
  "reagent",
  "catalyst",
  "substrate",
  "product",
  "reactant",
  "equilibrium",
  "thermodynamics",
  "kinetics",
  "enthalpy",
  "entropy",
  "stereochemistry",
  "isomer",
  "conformation",
  "resonance",
  "aromatic",
  "aliphatic",
  "heterocycle",
  "polymer",
  "biochemistry",
  "pharmaceutical",
  "drug",
  "medicinal",
];

/**
 * Check if a prompt is chemistry-related
 */
export function isChemistryPrompt(input: string): boolean {
  if (!input || input.trim().length === 0) {
    return false;
  }
  
  const inputLower = input.toLowerCase();
  return CHEMISTRY_KEYWORDS.some((keyword) => inputLower.includes(keyword));
}

/**
 * Get validation error message
 */
export function getValidationError(): string {
  return "I can help only with chemistry-related topics.";
}

