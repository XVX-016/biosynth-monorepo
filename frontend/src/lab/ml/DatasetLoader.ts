/**
 * DatasetLoader - Load and normalize molecular datasets
 */

import type { MoleculeState } from '../engines/MoleculeStateEngine';

export interface DatasetItem {
  id: string;
  state: MoleculeState;
  target: number | number[]; // Single value or multi-task
  metadata?: Record<string, any>;
}

export interface DatasetStats {
  count: number;
  mean?: number;
  std?: number;
  min?: number;
  max?: number;
}

/**
 * DatasetLoader class
 */
export class DatasetLoader {
  private normalizationParams: DatasetStats | null = null;

  /**
   * Initialize loader
   */
  async init(): Promise<void> {
    // Initialize any required resources
    console.log('DatasetLoader initialized');
  }

  /**
   * Load dataset from path/URL
   */
  async load(path: string): Promise<DatasetItem[]> {
    // TODO: Implement actual dataset loading
    // This would typically:
    // 1. Load from file (JSON, CSV, SDF, etc.)
    // 2. Parse molecule data
    // 3. Extract targets (properties, labels, etc.)
    // 4. Convert to DatasetItem format

    try {
      // Example: Load from JSON file
      const response = await fetch(path);
      if (!response.ok) {
        throw new Error(`Failed to load dataset: ${response.statusText}`);
      }

      const data = await response.json();
      
      // Convert to DatasetItem format
      const items: DatasetItem[] = data.map((item: any, idx: number) => ({
        id: item.id || `item-${idx}`,
        state: this.parseMoleculeState(item.molecule || item.state),
        target: item.target || item.value || 0,
        metadata: item.metadata || {},
      }));

      return items;
    } catch (error) {
      console.error('Failed to load dataset:', error);
      // Return empty dataset on error
      return [];
    }
  }

  /**
   * Parse molecule state from various formats
   */
  private parseMoleculeState(data: any): MoleculeState {
    // Convert from JSON format to MoleculeState
    const atoms = new Map();
    const bonds = new Map();

    if (data.atoms) {
      data.atoms.forEach((atom: any) => {
        atoms.set(atom.id, {
          id: atom.id,
          element: atom.element,
          x: atom.x || 0,
          y: atom.y || 0,
          z: atom.z || 0,
          charge: atom.charge || 0,
        });
      });
    }

    if (data.bonds) {
      data.bonds.forEach((bond: any) => {
        bonds.set(bond.id, {
          id: bond.id,
          atoms: bond.atoms,
          order: bond.order || 1,
        });
      });
    }

    return { atoms, bonds };
  }

  /**
   * Normalize dataset (standardize targets)
   */
  normalize(data: DatasetItem[]): DatasetItem[] {
    if (data.length === 0) return data;

    // Calculate statistics
    const targets = data.map((item) => {
      const target = Array.isArray(item.target) ? item.target[0] : item.target;
      return typeof target === 'number' ? target : 0;
    });

    const mean = targets.reduce((sum, t) => sum + t, 0) / targets.length;
    const variance =
      targets.reduce((sum, t) => sum + Math.pow(t - mean, 2), 0) / targets.length;
    const std = Math.sqrt(variance);
    const min = Math.min(...targets);
    const max = Math.max(...targets);

    this.normalizationParams = {
      count: data.length,
      mean,
      std,
      min,
      max,
    };

    // Normalize targets
    const normalized = data.map((item) => {
      const target = Array.isArray(item.target) ? item.target[0] : item.target;
      const numTarget = typeof target === 'number' ? target : 0;

      // Z-score normalization
      const normalizedTarget = std > 0 ? (numTarget - mean) / std : numTarget;

      return {
        ...item,
        target: Array.isArray(item.target)
          ? [normalizedTarget, ...item.target.slice(1)]
          : normalizedTarget,
      };
    });

    return normalized;
  }

  /**
   * Denormalize predictions
   */
  denormalize(predictions: number[]): number[] {
    if (!this.normalizationParams || !this.normalizationParams.mean || !this.normalizationParams.std) {
      return predictions;
    }

    const { mean, std } = this.normalizationParams;
    return predictions.map((pred) => pred * std + mean);
  }

  /**
   * Get dataset statistics
   */
  getStats(): DatasetStats | null {
    return this.normalizationParams;
  }

  /**
   * Split dataset into train/validation/test
   */
  split(
    data: DatasetItem[],
    trainRatio: number = 0.7,
    valRatio: number = 0.15
  ): {
    train: DatasetItem[];
    validation: DatasetItem[];
    test: DatasetItem[];
  } {
    const shuffled = [...data].sort(() => Math.random() - 0.5);
    const trainEnd = Math.floor(shuffled.length * trainRatio);
    const valEnd = trainEnd + Math.floor(shuffled.length * valRatio);

    return {
      train: shuffled.slice(0, trainEnd),
      validation: shuffled.slice(trainEnd, valEnd),
      test: shuffled.slice(valEnd),
    };
  }
}

