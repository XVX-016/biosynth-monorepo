/**
 * TrainLoop - Training loop with checkpointing
 */

import type { GNNModel } from './GNNModel';
import type { GraphFeature } from './Featurizer';
import type { TrainingConfig, TrainingMetrics } from './TrainingPipeline';

export interface TrainingData {
  graph: GraphFeature;
  target: number | number[];
}

export interface Checkpoint {
  epoch: number;
  model: GNNModel;
  metrics: TrainingMetrics;
  timestamp: number;
}

/**
 * TrainLoop class
 */
export class TrainLoop {
  private checkpoints: Checkpoint[] = [];
  private bestLoss: number = Infinity;
  private bestEpoch: number = 0;

  /**
   * Train model
   */
  async train(
    model: GNNModel,
    trainData: TrainingData[],
    valData: TrainingData[],
    config: TrainingConfig
  ): Promise<{
    model: GNNModel;
    metrics: TrainingMetrics[];
    bestEpoch: number;
    finalLoss: number;
  }> {
    const metrics: TrainingMetrics[] = [];

    // Initialize model if needed
    if (!model.getConfig()) {
      await model.init();
    }

    // Training loop
    for (let epoch = 0; epoch < config.epochs; epoch++) {
      // Shuffle training data
      const shuffled = [...trainData].sort(() => Math.random() - 0.5);

      // Batch training
      const batchSize = config.batchSize;
      let epochLoss = 0;
      let batchCount = 0;

      for (let i = 0; i < shuffled.length; i += batchSize) {
        const batch = shuffled.slice(i, i + batchSize);
        const batchLoss = await this.trainBatch(model, batch, config.learningRate);
        epochLoss += batchLoss;
        batchCount++;
      }

      const avgLoss = epochLoss / batchCount;

      // Validation
      let valLoss = 0;
      if (valData.length > 0) {
        const valPredictions = await model.forwardBatch(
          valData.map((d) => d.graph)
        );
        valLoss = this.computeLoss(
          valPredictions,
          valData.map((d) => (Array.isArray(d.target) ? d.target[0] : d.target))
        );
      }

      // Record metrics
      const metric: TrainingMetrics = {
        epoch,
        loss: avgLoss,
        validationLoss: valData.length > 0 ? valLoss : undefined,
      };

      metrics.push(metric);

      // Check for best model
      const currentLoss = valData.length > 0 ? valLoss : avgLoss;
      if (currentLoss < this.bestLoss) {
        this.bestLoss = currentLoss;
        this.bestEpoch = epoch;
      }

      // Save checkpoint
      if (epoch % config.checkpointInterval === 0 || epoch === config.epochs - 1) {
        await this.saveCheckpoint(model, metric, epoch);
      }

      console.log(
        `Epoch ${epoch}/${config.epochs}: Loss=${avgLoss.toFixed(4)}, ValLoss=${valLoss.toFixed(4)}`
      );
    }

    return {
      model,
      metrics,
      bestEpoch: this.bestEpoch,
      finalLoss: metrics[metrics.length - 1].loss,
    };
  }

  /**
   * Train on a single batch
   */
  private async trainBatch(
    model: GNNModel,
    batch: TrainingData[],
    learningRate: number
  ): Promise<number> {
    // Forward pass
    const predictions = await model.forwardBatch(batch.map((d) => d.graph));
    const targets = batch.map((d) => (Array.isArray(d.target) ? d.target[0] : d.target));

    // Compute loss
    const loss = this.computeLoss(predictions, targets);

    // TODO: Backward pass and weight update
    // This requires actual gradient computation and optimizer
    // For now, this is a placeholder

    return loss;
  }

  /**
   * Compute loss (MSE for regression)
   */
  private computeLoss(predictions: number[], targets: number[]): number {
    if (predictions.length !== targets.length) {
      return Infinity;
    }

    let sumSquaredError = 0;
    for (let i = 0; i < predictions.length; i++) {
      const error = predictions[i] - targets[i];
      sumSquaredError += error * error;
    }

    return sumSquaredError / predictions.length;
  }

  /**
   * Save checkpoint
   */
  async saveCheckpoint(
    model: GNNModel,
    metrics: TrainingMetrics,
    epoch: number
  ): Promise<void> {
    const checkpoint: Checkpoint = {
      epoch,
      model,
      metrics,
      timestamp: Date.now(),
    };

    this.checkpoints.push(checkpoint);

    // Save to storage (localStorage or backend)
    try {
      const checkpointData = {
        epoch,
        metrics,
        timestamp: checkpoint.timestamp,
        // Model weights would be serialized here
      };

      localStorage.setItem(`model-checkpoint-${epoch}`, JSON.stringify(checkpointData));
    } catch (error) {
      console.warn('Failed to save checkpoint to localStorage:', error);
    }

    console.log(`Checkpoint saved at epoch ${epoch}`);
  }

  /**
   * Load checkpoint
   */
  async loadCheckpoint(epoch: number): Promise<Checkpoint | null> {
    // Try to load from checkpoints array first
    const checkpoint = this.checkpoints.find((c) => c.epoch === epoch);
    if (checkpoint) {
      return checkpoint;
    }

    // Try to load from storage
    try {
      const data = localStorage.getItem(`model-checkpoint-${epoch}`);
      if (data) {
        const checkpointData = JSON.parse(data);
        // Would reconstruct model from weights here
        return null; // Placeholder
      }
    } catch (error) {
      console.warn('Failed to load checkpoint:', error);
    }

    return null;
  }

  /**
   * Get all checkpoints
   */
  getCheckpoints(): Checkpoint[] {
    return [...this.checkpoints];
  }

  /**
   * Get best checkpoint
   */
  getBestCheckpoint(): Checkpoint | null {
    if (this.checkpoints.length === 0) return null;

    return this.checkpoints.find((c) => c.epoch === this.bestEpoch) || this.checkpoints[0];
  }
}

export const trainLoop = new TrainLoop();

