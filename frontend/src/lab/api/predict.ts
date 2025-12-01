/**
 * Predict API - Backend endpoint for model predictions
 * 
 * Note: This is a frontend API client
 * The actual endpoint should be implemented in backend
 */

import type { MoleculeState } from '../engines/MoleculeStateEngine';

export interface PredictionRequest {
  molecule?: MoleculeState;
  smiles?: string;
  graph?: any;
  node_order?: string[];
  modelId?: string;
  properties?: string[]; // e.g., ['logP', 'toxicity', 'solubility']
  return_attention?: boolean;
  attention_layer?: string; // 'last', 'all', or layer index
}

export interface PredictionResponse {
  predictions: Record<string, number>;
  model_id: string;
  confidence?: Record<string, number>;
  attentions?: {
    [key: string]: number[];
    node_importance?: number[];
  };
  edge_index?: number[][];
  node_mapping?: Record<number, string>;
  warnings?: string[];
}

/**
 * Convert MoleculeState to graph JSON format
 */
function moleculeToGraph(molecule: MoleculeState): any {
  const atoms = Array.from(molecule.atoms.values()).map((atom) => ({
    id: atom.id,
    element: atom.element,
    x: atom.x,
    y: atom.y,
    z: 0, // Default z if not available
    charge: atom.charge || 0,
  }));

  const bonds = Array.from(molecule.bonds.values()).map((bond) => ({
    id: bond.id,
    atoms: bond.atoms,
    order: bond.order,
  }));

  const node_order = atoms.map((a) => a.id);

  return { atoms, bonds, node_order };
}

/**
 * Predict molecule properties
 */
export async function predict(
  request: PredictionRequest
): Promise<PredictionResponse> {
  try {
    // Build request payload
    const payload: any = {
      model_id: request.modelId,
      properties: request.properties || ['logP', 'toxicity', 'solubility'],
      return_attention: request.return_attention || false,
      attention_layer: request.attention_layer || 'last',
    };

    if (request.smiles) {
      payload.smiles = request.smiles;
      if (request.node_order) {
        payload.node_order = request.node_order;
      }
    } else if (request.graph) {
      payload.graph = request.graph;
    } else if (request.molecule) {
      payload.graph = moleculeToGraph(request.molecule);
    } else {
      throw new Error('No molecule input provided');
    }

    // Send to backend
    const response = await fetch('/api/predict/property', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify(payload),
    });

    if (!response.ok) {
      const errorData = await response.json().catch(() => ({}));
      throw new Error(errorData.detail || `Prediction failed: ${response.statusText}`);
    }

    const data = await response.json();
    return data;
  } catch (error) {
    console.error('Prediction error:', error);
    throw error;
  }
}

/**
 * Batch prediction
 */
export async function predictBatch(
  molecules: MoleculeState[],
  modelId?: string
): Promise<PredictionResponse[]> {
  const featurizer = new Featurizer();
  const features = molecules.map((mol) => ({
    graph: featurizer.computeGraph(mol),
    ecfp: featurizer.computeECFP(mol).fingerprint,
  }));

  try {
    const response = await fetch('/api/predict/batch', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        features,
        modelId: modelId || 'default',
      }),
    });

    if (!response.ok) {
      throw new Error(`Batch prediction failed: ${response.statusText}`);
    }

    const data = await response.json();
    return data.predictions;
  } catch (error) {
    console.error('Batch prediction error:', error);
    throw error;
  }
}

