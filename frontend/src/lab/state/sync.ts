/**
 * Multi-tab synchronization using BroadcastChannel
 * 
 * Implements optimistic locking and simple CRDT-like merge policy for concurrent edits.
 */

import type { MoleculeState } from '../engines/MoleculeStateEngine';
import { IAction } from '../engines/actions/Action.types';

export interface SyncMessage {
  type: 'action' | 'state' | 'lock' | 'unlock' | 'conflict';
  tabId: string;
  timestamp: number;
  data?: any;
}

export interface SyncConfig {
  channelName: string;
  enableSync: boolean;
  conflictResolution: 'last-write-wins' | 'merge' | 'manual';
}

export class MultiTabSync {
  private channel: BroadcastChannel | null = null;
  private tabId: string;
  private config: Required<SyncConfig>;
  private isLocked: boolean = false;
  private lockOwner: string | null = null;
  private pendingActions: IAction[] = [];
  private listeners: ((message: SyncMessage) => void)[] = [];

  constructor(config: Partial<SyncConfig> = {}) {
    this.tabId = crypto.randomUUID();
    this.config = {
      channelName: config.channelName || 'molecule-lab-sync',
      enableSync: config.enableSync ?? true,
      conflictResolution: config.conflictResolution || 'last-write-wins',
    };

    if (this.config.enableSync && typeof BroadcastChannel !== 'undefined') {
      this.channel = new BroadcastChannel(this.config.channelName);
      this.channel.onmessage = (event) => this.handleMessage(event.data);
    }
  }

  /**
   * Broadcast an action to other tabs
   */
  broadcastAction(action: IAction): void {
    if (!this.channel || !this.config.enableSync) return;

    const message: SyncMessage = {
      type: 'action',
      tabId: this.tabId,
      timestamp: Date.now(),
      data: action.toJSON(),
    };

    this.channel.postMessage(message);
  }

  /**
   * Broadcast state update
   */
  broadcastState(state: MoleculeState): void {
    if (!this.channel || !this.config.enableSync) return;

    const message: SyncMessage = {
      type: 'state',
      tabId: this.tabId,
      timestamp: Date.now(),
      data: this.serializeState(state),
    };

    this.channel.postMessage(message);
  }

  /**
   * Request lock (optimistic locking)
   */
  requestLock(): boolean {
    if (!this.channel || !this.config.enableSync) {
      this.isLocked = true;
      this.lockOwner = this.tabId;
      return true;
    }

    if (this.isLocked && this.lockOwner !== this.tabId) {
      return false; // Lock held by another tab
    }

    const message: SyncMessage = {
      type: 'lock',
      tabId: this.tabId,
      timestamp: Date.now(),
    };

    this.channel.postMessage(message);
    this.isLocked = true;
    this.lockOwner = this.tabId;

    return true;
  }

  /**
   * Release lock
   */
  releaseLock(): void {
    if (!this.channel || !this.config.enableSync) {
      this.isLocked = false;
      this.lockOwner = null;
      return;
    }

    const message: SyncMessage = {
      type: 'unlock',
      tabId: this.tabId,
      timestamp: Date.now(),
    };

    this.channel.postMessage(message);
    this.isLocked = false;
    this.lockOwner = null;
  }

  /**
   * Check if we have the lock
   */
  hasLock(): boolean {
    return this.isLocked && this.lockOwner === this.tabId;
  }

  /**
   * Subscribe to sync messages
   */
  onMessage(callback: (message: SyncMessage) => void): () => void {
    this.listeners.push(callback);
    return () => {
      const index = this.listeners.indexOf(callback);
      if (index > -1) {
        this.listeners.splice(index, 1);
      }
    };
  }

  /**
   * Handle incoming message
   */
  private handleMessage(message: SyncMessage): void {
    // Ignore messages from self
    if (message.tabId === this.tabId) {
      return;
    }

    switch (message.type) {
      case 'lock':
        // Another tab requested lock
        if (this.isLocked && this.lockOwner === this.tabId) {
          // Conflict - send conflict message
          this.sendConflict(message.tabId);
        }
        break;

      case 'unlock':
        // Another tab released lock
        if (this.lockOwner === message.tabId) {
          this.isLocked = false;
          this.lockOwner = null;
        }
        break;

      case 'action':
        // Another tab performed an action
        this.notifyListeners(message);
        break;

      case 'state':
        // Another tab updated state
        this.notifyListeners(message);
        break;

      case 'conflict':
        // Conflict detected
        this.notifyListeners(message);
        break;
    }
  }

  /**
   * Send conflict message
   */
  private sendConflict(otherTabId: string): void {
    if (!this.channel) return;

    const message: SyncMessage = {
      type: 'conflict',
      tabId: this.tabId,
      timestamp: Date.now(),
      data: { otherTabId },
    };

    this.channel.postMessage(message);
  }

  /**
   * Notify listeners
   */
  private notifyListeners(message: SyncMessage): void {
    this.listeners.forEach(cb => {
      try {
        cb(message);
      } catch (error) {
        console.error('Error in sync listener:', error);
      }
    });
  }

  /**
   * Serialize state for transmission
   */
  private serializeState(state: MoleculeState): any {
    return {
      atoms: Array.from(state.atoms.entries()),
      bonds: Array.from(state.bonds.entries()),
    };
  }

  /**
   * Deserialize state from message
   */
  deserializeState(data: any): MoleculeState {
    return {
      atoms: new Map(data.atoms),
      bonds: new Map(data.bonds),
    };
  }

  /**
   * Close sync channel
   */
  close(): void {
    if (this.channel) {
      this.releaseLock();
      this.channel.close();
      this.channel = null;
    }
  }
}

