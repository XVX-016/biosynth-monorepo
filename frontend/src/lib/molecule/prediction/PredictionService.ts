/**
 * PredictionService.ts
 *
 * Phase 7+: ML Prediction Pipeline
 *
 * Provides debounced, validated prediction requests with:
 * - Normalized SMILES input
 * - Validation integration
 * - Async queue (300ms debounce)
 * - Attention map support
 * - Error handling
 */

import type { Molecule } from '../Molecule'
import { toSMILES } from '../export'
import { validateMolecule } from '../validation/Validator'

export interface AttentionData {
  edgeAttentions: number[]
  edgeIndex: number[][]
  nodeImportance?: number[]
}

export interface PredictionResult {
  properties: Record<string, number>
  attention?: AttentionData
  modelId?: string
  timestamp: number
}

export interface PredictionRequest {
  molecule: Molecule
  properties?: string[]
  modelId?: string
  returnAttention?: boolean
}

export class PredictionService {
  private pendingRequest: NodeJS.Timeout | null = null
  private lastMoleculeHash: string | null = null
  private requestQueue: PredictionRequest[] = []
  private isProcessing = false
  private readonly DEBOUNCE_MS = 300

  /**
   * Predict properties for a molecule (debounced)
   */
  async predict(request: PredictionRequest): Promise<PredictionResult> {
    const validation = validateMolecule(request.molecule)
    if (!validation.valid && validation.errors.length > 0) {
      throw new Error(`Invalid molecule: ${validation.errors[0].message}`)
    }

    if (this.pendingRequest) {
      clearTimeout(this.pendingRequest)
      this.pendingRequest = null
    }

    this.requestQueue.push(request)

    return new Promise((resolve, reject) => {
      this.pendingRequest = setTimeout(async () => {
        try {
          const result = await this.processQueue()
          resolve(result)
        } catch (error) {
          reject(error)
        }
      }, this.DEBOUNCE_MS)
    })
  }

  /**
   * Process the prediction queue
   */
  private async processQueue(): Promise<PredictionResult> {
    if (this.isProcessing || this.requestQueue.length === 0) {
      throw new Error('No requests in queue')
    }

    this.isProcessing = true

    try {
      const request = this.requestQueue[this.requestQueue.length - 1]
      this.requestQueue = []

      const smiles = await toSMILES(request.molecule, true)
      if (!smiles) throw new Error('Failed to generate SMILES')

      const moleculeHash = this.hashMolecule(request.molecule)
      if (moleculeHash === this.lastMoleculeHash && this.lastMoleculeHash !== null) {
        // For now, still process request
      }
      this.lastMoleculeHash = moleculeHash

      const requestBody = {
        smiles,
        properties: request.properties,
        model_id: request.modelId,
        return_attention: request.returnAttention || false,
      }

      const response = await fetch('/api/predict/property', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(requestBody),
      })

      if (!response.ok) {
        const errorData = await response.json().catch(() => ({ detail: 'Unknown error' }))
        throw new Error(errorData.detail || `Prediction failed: ${response.statusText}`)
      }

      const data = await response.json()

      const result: PredictionResult = {
        properties: data.predictions || {},
        modelId: data.model_id,
        timestamp: Date.now(),
      }

      if (request.returnAttention && data.attention) {
        result.attention = {
          edgeAttentions: data.attention.edge_attentions || [],
          edgeIndex: data.attention.edge_index || [],
          nodeImportance: data.attention.node_importance,
        }
      }

      return result
    } finally {
      this.isProcessing = false
    }
  }

  /**
   * Predict with attention map
   */
  async predictWithAttention(
    molecule: Molecule,
    modelId?: string,
    properties?: string[]
  ): Promise<PredictionResult> {
    return this.predict({
      molecule,
      properties,
      modelId,
      returnAttention: true,
    })
  }

  /**
   * Batch prediction (no debounce)
   */
  async predictBatch(
    molecules: Molecule[],
    modelId?: string,
    properties?: string[]
  ): Promise<PredictionResult[]> {
    for (const mol of molecules) {
      const validation = validateMolecule(mol)
      if (!validation.valid && validation.errors.length > 0) {
        throw new Error(`Invalid molecule: ${validation.errors[0].message}`)
      }
    }

    const inputs = await Promise.all(
      molecules.map(async (mol) => {
        const smiles = await toSMILES(mol, true)
        if (!smiles) throw new Error('Failed to generate SMILES')
        return { smiles }
      })
    )

    const response = await fetch('/api/predict/batch', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ inputs, model_id: modelId, properties }),
    })

    if (!response.ok) {
      const errorData = await response.json().catch(() => ({ detail: 'Unknown error' }))
      throw new Error(errorData.detail || `Batch prediction failed: ${response.statusText}`)
    }

    const data = await response.json()

    return data.results.map((r: any) => ({
      properties: r.predictions || {},
      modelId: r.model_id,
      timestamp: Date.now(),
    }))
  }

  /**
   * Cancel pending predictions
   */
  cancel(): void {
    if (this.pendingRequest) clearTimeout(this.pendingRequest)
    this.pendingRequest = null
    this.requestQueue = []
  }

  /**
   * Hash molecule for caching
   */
  private hashMolecule(molecule: Molecule): string {
    const atoms = molecule
      .getAtoms()
      .map((a) => `${a.element}:${a.position.join(',')}`)
      .sort()
      .join('|')
    const bonds = molecule
      .getBonds()
      .map((b) => `${b.atom1}-${b.atom2}:${b.order}`)
      .sort()
      .join('|')
    return `${atoms}|${bonds}`
  }
}

// Singleton instance
export const predictionService = new PredictionService()
