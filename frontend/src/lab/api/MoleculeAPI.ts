/**
 * MoleculeAPI - Backend API integration for molecule operations
 */

import type { MoleculeState } from '../engines/MoleculeStateEngine';

export interface SMILESResponse {
  smiles: string;
  success: boolean;
  error?: string;
}

export interface SDFResponse {
  sdf: string;
  success: boolean;
  error?: string;
}

export interface ValidationResponse {
  valid: boolean;
  issues: Array<{
    type: string;
    message: string;
    atomId?: string;
    bondId?: string;
  }>;
  sanitized?: MoleculeState;
}

/**
 * Convert molecule state to SMILES via backend
 */
export async function convertToSMILES(state: MoleculeState): Promise<SMILESResponse> {
  try {
    // Convert state to format backend expects
    const payload = {
      atoms: Array.from(state.atoms.values()).map((atom) => ({
        id: atom.id,
        element: atom.element,
        x: atom.x,
        y: atom.y,
        z: atom.z || 0,
        charge: atom.charge,
      })),
      bonds: Array.from(state.bonds.values()).map((bond) => ({
        id: bond.id,
        atoms: bond.atoms,
        order: bond.order,
      })),
    };

    const response = await fetch('/api/mol/to-smiles', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify(payload),
    });

    if (!response.ok) {
      throw new Error(`HTTP error! status: ${response.status}`);
    }

    const data = await response.json();
    return {
      smiles: data.smiles || '',
      success: data.success !== false,
      error: data.error,
    };
  } catch (error) {
    console.error('Failed to convert to SMILES:', error);
    return {
      smiles: '',
      success: false,
      error: error instanceof Error ? error.message : 'Unknown error',
    };
  }
}

/**
 * Convert molecule state to SDF via backend
 */
export async function convertToSDF(state: MoleculeState): Promise<SDFResponse> {
  try {
    const payload = {
      atoms: Array.from(state.atoms.values()).map((atom) => ({
        id: atom.id,
        element: atom.element,
        x: atom.x,
        y: atom.y,
        z: atom.z || 0,
        charge: atom.charge,
      })),
      bonds: Array.from(state.bonds.values()).map((bond) => ({
        id: bond.id,
        atoms: bond.atoms,
        order: bond.order,
      })),
    };

    const response = await fetch('/api/mol/to-sdf', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify(payload),
    });

    if (!response.ok) {
      throw new Error(`HTTP error! status: ${response.status}`);
    }

    const data = await response.json();
    return {
      sdf: data.sdf || '',
      success: data.success !== false,
      error: data.error,
    };
  } catch (error) {
    console.error('Failed to convert to SDF:', error);
    return {
      sdf: '',
      success: false,
      error: error instanceof Error ? error.message : 'Unknown error',
    };
  }
}

/**
 * Validate molecule via backend (deep validation with RDKit)
 */
export async function validateMolecule(state: MoleculeState): Promise<ValidationResponse> {
  try {
    const payload = {
      atoms: Array.from(state.atoms.values()),
      bonds: Array.from(state.bonds.values()),
    };

    const response = await fetch('/api/mol/validate', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify(payload),
    });

    if (!response.ok) {
      throw new Error(`HTTP error! status: ${response.status}`);
    }

    const data = await response.json();
    return {
      valid: data.valid !== false,
      issues: data.issues || [],
      sanitized: data.sanitized,
    };
  } catch (error) {
    console.error('Failed to validate molecule:', error);
    return {
      valid: false,
      issues: [
        {
          type: 'network_error',
          message: error instanceof Error ? error.message : 'Validation failed',
        },
      ],
    };
  }
}

/**
 * Parse SMILES and return molecule state
 */
export async function parseSMILES(smiles: string): Promise<MoleculeState | null> {
  try {
    const response = await fetch('/api/mol/parse-smiles', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ smiles }),
    });

    if (!response.ok) {
      throw new Error(`HTTP error! status: ${response.status}`);
    }

    const data = await response.json();
    if (!data.success || !data.molecule) {
      return null;
    }

    // Convert backend format to MoleculeState
    const atoms = new Map();
    const bonds = new Map();

    data.molecule.atoms?.forEach((atom: any) => {
      atoms.set(atom.id, {
        id: atom.id,
        element: atom.element,
        x: atom.x || 0,
        y: atom.y || 0,
        z: atom.z || 0,
        charge: atom.charge || 0,
      });
    });

    data.molecule.bonds?.forEach((bond: any) => {
      bonds.set(bond.id, {
        id: bond.id,
        atoms: bond.atoms,
        order: bond.order || 1,
      });
    });

    return { atoms, bonds };
  } catch (error) {
    console.error('Failed to parse SMILES:', error);
    return null;
  }
}

