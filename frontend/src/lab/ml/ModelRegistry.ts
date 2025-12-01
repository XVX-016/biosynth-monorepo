/**
 * ModelRegistry - Registry for ML models
 */

import type { GNNModel } from './GNNModel';

export interface ModelInfo {
  id: string;
  name: string;
  description?: string;
  type: 'gnn' | 'ecfp' | 'hybrid';
  version: string;
  path?: string;
  metadata?: Record<string, any>;
  createdAt: number;
  updatedAt: number;
}

export interface ModelMetrics {
  accuracy?: number;
  loss?: number;
  validationAccuracy?: number;
  validationLoss?: number;
  trainingTime?: number;
}

/**
 * ModelRegistry class
 */
export class ModelRegistry {
  private models = new Map<string, ModelInfo>();
  private loadedModels = new Map<string, GNNModel>();

  /**
   * Register a model
   */
  register(model: ModelInfo): void {
    this.models.set(model.id, {
      ...model,
      updatedAt: Date.now(),
    });
  }

  /**
   * Get model info
   */
  get(id: string): ModelInfo | undefined {
    return this.models.get(id);
  }

  /**
   * List all models
   */
  listModels(): ModelInfo[] {
    return Array.from(this.models.values());
  }

  /**
   * List models by type
   */
  listByType(type: ModelInfo['type']): ModelInfo[] {
    return Array.from(this.models.values()).filter((m) => m.type === type);
  }

  /**
   * Load model into memory
   */
  async loadModel(id: string): Promise<GNNModel | null> {
    // Check if already loaded
    if (this.loadedModels.has(id)) {
      return this.loadedModels.get(id)!;
    }

    const modelInfo = this.models.get(id);
    if (!modelInfo) {
      console.error(`Model ${id} not found in registry`);
      return null;
    }

    try {
      // Load model from path or backend
      if (modelInfo.path) {
        const model = await GNNModel.load(modelInfo.path);
        this.loadedModels.set(id, model);
        return model;
      } else {
        // Load from backend API
        const response = await fetch(`/api/models/${id}/load`);
        if (!response.ok) {
          throw new Error(`Failed to load model: ${response.statusText}`);
        }
        // Model would be reconstructed from response
        // For now, return null
        return null;
      }
    } catch (error) {
      console.error(`Failed to load model ${id}:`, error);
      return null;
    }
  }

  /**
   * Unload model from memory
   */
  unloadModel(id: string): void {
    this.loadedModels.delete(id);
  }

  /**
   * Remove model from registry
   */
  remove(id: string): boolean {
    this.unloadModel(id);
    return this.models.delete(id);
  }

  /**
   * Update model info
   */
  update(id: string, updates: Partial<ModelInfo>): boolean {
    const model = this.models.get(id);
    if (!model) {
      return false;
    }

    this.models.set(id, {
      ...model,
      ...updates,
      updatedAt: Date.now(),
    });

    return true;
  }

  /**
   * Get default model
   */
  getDefault(): ModelInfo | null {
    const models = this.listModels();
    if (models.length === 0) return null;

    // Try to find a model marked as default
    const defaultModel = models.find((m) => m.metadata?.default === true);
    if (defaultModel) return defaultModel;

    // Return most recently updated
    return models.sort((a, b) => b.updatedAt - a.updatedAt)[0];
  }

  /**
   * Initialize with default models
   */
  async initDefaults(): Promise<void> {
    // Register default models
    this.register({
      id: 'default-gnn',
      name: 'Default GNN',
      description: 'Default graph neural network for property prediction',
      type: 'gnn',
      version: '1.0.0',
      createdAt: Date.now(),
      updatedAt: Date.now(),
    });

    this.register({
      id: 'default-ecfp',
      name: 'Default ECFP',
      description: 'ECFP-based property predictor',
      type: 'ecfp',
      version: '1.0.0',
      createdAt: Date.now(),
      updatedAt: Date.now(),
    });
  }
}

export const modelRegistry = new ModelRegistry();

