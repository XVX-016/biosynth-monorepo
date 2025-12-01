/**
 * predictWithAttention - API client for attention-based predictions
 */

import type { MoleculeState } from '../engines/MoleculeStateEngine';
import { predict, type PredictionRequest, type PredictionResponse } from './predict';

export interface AttentionPredictRequest {
  molecule?: MoleculeState;
  smiles?: string;
  graph?: any;
  node_order?: string[];
  modelId?: string;
  returnAttention?: boolean;
  attention_layer?: string; // 'last', 'all', or layer index
}

export interface AttentionResult {
  prediction: any;
  attentions: number[][];
  edgeIndex: number[][];
  nodeImportance?: number[];
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
    z: 0,
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
 * Predict with attention weights
 */
export async function predictWithAttention(
  request: AttentionPredictRequest
): Promise<AttentionResult> {
  try {
    const payload: PredictionRequest = {
      molecule: request.molecule,
      smiles: request.smiles,
      graph: request.graph || (request.molecule ? moleculeToGraph(request.molecule) : undefined),
      node_order: request.node_order,
      modelId: request.modelId,
      return_attention: request.returnAttention !== false,
      attention_layer: request.attention_layer || 'last',
    };

    const data = await predict(payload);

    // Extract attention data
    const attentions = data.attentions || {};
    const edgeIndex = data.edge_index || [];
    
    // Get the requested layer (default to 'layer_last')
    const layerKey = request.attention_layer === 'all' 
      ? Object.keys(attentions)[0] 
      : `layer_${request.attention_layer || 'last'}`;
    
    const edgeAttentions = attentions[layerKey] || attentions['layer_last'] || [];
    const nodeImportance = attentions['node_importance'];

    return {
      prediction: data.predictions,
      attentions: request.attention_layer === 'all' 
        ? Object.values(attentions).filter((v): v is number[] => Array.isArray(v))
        : [edgeAttentions],
      edgeIndex: edgeIndex,
      nodeImportance: nodeImportance,
    };
  } catch (error) {
    console.error('Attention prediction error:', error);
    throw error;
  }
}

/**
 * Batch prediction with attention
 */
export async function predictBatchWithAttention(
  molecules: MoleculeState[],
  modelId?: string
): Promise<AttentionResult[]> {
  const featurizer = new Featurizer();
  const features = molecules.map((mol) => ({
    graph: featurizer.computeGraph(mol),
  }));

  try {
    const response = await fetch('/api/ml/predict/batch', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        features,
        modelId: modelId || 'attention-gnn',
        returnAttention: true,
      }),
    });

    if (!response.ok) {
      throw new Error(`Batch prediction failed: ${response.statusText}`);
    }

    const data = await response.json();
    return data.predictions;
  } catch (error) {
    console.error('Batch attention prediction error:', error);
    throw error;
  }
}

