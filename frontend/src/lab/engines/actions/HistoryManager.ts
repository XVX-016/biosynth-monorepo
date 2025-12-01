/**
 * HistoryManager - Advanced history management with compression and snapshots
 * 
 * Manages undo/redo history with size limits, compression, and snapshot isolation.
 */

import { IAction } from './Action.types';

export interface HistorySnapshot {
  id: string;
  timestamp: number;
  label: string;
  state: any;
  actionCount: number;
}

export interface HistoryConfig {
  maxSize: number;
  compressionThreshold: number;
  snapshotInterval: number; // Create snapshot every N actions
  enableCompression: boolean;
}

export class HistoryManager {
  private undoStack: IAction[] = [];
  private redoStack: IAction[] = [];
  private snapshots: HistorySnapshot[] = [];
  private config: Required<HistoryConfig>;
  private actionCountSinceSnapshot = 0;

  constructor(config: Partial<HistoryConfig> = {}) {
    this.config = {
      maxSize: config.maxSize ?? 100,
      compressionThreshold: config.compressionThreshold ?? 50,
      snapshotInterval: config.snapshotInterval ?? 20,
      enableCompression: config.enableCompression ?? true,
    };
  }

  /**
   * Add action to history
   */
  push(action: IAction, state: any): void {
    this.undoStack.push(action);
    this.redoStack = []; // Clear redo stack on new action

    // Check if we need to create a snapshot
    this.actionCountSinceSnapshot++;
    if (this.actionCountSinceSnapshot >= this.config.snapshotInterval) {
      this.createSnapshot(state, `Snapshot after ${this.undoStack.length} actions`);
      this.actionCountSinceSnapshot = 0;
    }

    // Trim history if needed
    this.trimHistory();
  }

  /**
   * Pop action for undo
   */
  pop(): IAction | undefined {
    const action = this.undoStack.pop();
    if (action) {
      this.actionCountSinceSnapshot--;
    }
    return action;
  }

  /**
   * Push to redo stack
   */
  pushRedo(action: IAction): void {
    this.redoStack.push(action);
  }

  /**
   * Pop from redo stack
   */
  popRedo(): IAction | undefined {
    return this.redoStack.pop();
  }

  /**
   * Create a snapshot of current state
   */
  createSnapshot(state: any, label?: string): HistorySnapshot {
    const snapshot: HistorySnapshot = {
      id: crypto.randomUUID(),
      timestamp: Date.now(),
      label: label || `Snapshot at ${new Date().toISOString()}`,
      state: this.deepClone(state),
      actionCount: this.undoStack.length,
    };

    this.snapshots.push(snapshot);

    // Limit snapshot count
    if (this.snapshots.length > 10) {
      this.snapshots.shift();
    }

    return snapshot;
  }

  /**
   * Get nearest snapshot before given action count
   */
  getNearestSnapshot(actionCount: number): HistorySnapshot | null {
    let nearest: HistorySnapshot | null = null;
    let minDiff = Infinity;

    for (const snapshot of this.snapshots) {
      const diff = actionCount - snapshot.actionCount;
      if (diff >= 0 && diff < minDiff) {
        minDiff = diff;
        nearest = snapshot;
      }
    }

    return nearest;
  }

  /**
   * Compress history by merging small actions
   */
  compress(): number {
    if (!this.config.enableCompression || this.undoStack.length < this.config.compressionThreshold) {
      return 0;
    }

    // Simple compression: keep first 10%, middle 10%, and last 50%
    const keepCount = Math.floor(this.undoStack.length * 0.1);
    const lastKeepCount = Math.floor(this.undoStack.length * 0.5);

    const compressed: IAction[] = [];

    // Keep first N
    compressed.push(...this.undoStack.slice(0, keepCount));

    // Keep middle N (skip some)
    const middleStart = Math.floor(this.undoStack.length * 0.4);
    compressed.push(...this.undoStack.slice(middleStart, middleStart + keepCount));

    // Keep last N
    compressed.push(...this.undoStack.slice(-lastKeepCount));

    const removed = this.undoStack.length - compressed.length;
    this.undoStack = compressed;

    return removed;
  }

  /**
   * Trim history to max size
   */
  private trimHistory(): void {
    if (this.undoStack.length > this.config.maxSize) {
      const excess = this.undoStack.length - this.config.maxSize;
      
      // Try compression first
      if (this.config.enableCompression) {
        this.compress();
      }

      // If still too large, remove oldest
      if (this.undoStack.length > this.config.maxSize) {
        const toRemove = this.undoStack.length - this.config.maxSize;
        this.undoStack.splice(0, toRemove);
      }
    }
  }

  /**
   * Clear all history
   */
  clear(): void {
    this.undoStack = [];
    this.redoStack = [];
    this.snapshots = [];
    this.actionCountSinceSnapshot = 0;
  }

  /**
   * Get history statistics
   */
  getStats(): {
    undoCount: number;
    redoCount: number;
    snapshotCount: number;
    totalActions: number;
  } {
    return {
      undoCount: this.undoStack.length,
      redoCount: this.redoStack.length,
      snapshotCount: this.snapshots.length,
      totalActions: this.undoStack.length + this.redoStack.length,
    };
  }

  /**
   * Get undo stack (read-only)
   */
  getUndoStack(): readonly IAction[] {
    return [...this.undoStack];
  }

  /**
   * Get redo stack (read-only)
   */
  getRedoStack(): readonly IAction[] {
    return [...this.redoStack];
  }

  /**
   * Get snapshots (read-only)
   */
  getSnapshots(): readonly HistorySnapshot[] {
    return [...this.snapshots];
  }

  /**
   * Deep clone state for snapshots
   */
  private deepClone(obj: any): any {
    if (obj === null || typeof obj !== 'object') {
      return obj;
    }

    if (obj instanceof Map) {
      return new Map(Array.from(obj.entries()).map(([k, v]) => [k, this.deepClone(v)]));
    }

    if (obj instanceof Set) {
      return new Set(Array.from(obj).map(v => this.deepClone(v)));
    }

    if (Array.isArray(obj)) {
      return obj.map(v => this.deepClone(v));
    }

    const cloned: any = {};
    for (const key in obj) {
      if (obj.hasOwnProperty(key)) {
        cloned[key] = this.deepClone(obj[key]);
      }
    }

    return cloned;
  }
}

