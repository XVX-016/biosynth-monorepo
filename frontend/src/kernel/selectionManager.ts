/**
 * selectionManager.ts
 * Simple selection manager for the kernel. Not DOM-dependent.
 * Emits events: 'select', 'deselect', 'change'
 * 
 * This is a kernel-level selection manager, separate from the UI SelectionManager
 * in components/r3f/SelectionManager.ts which handles hover/drag state.
 */

type Handler = (payload?: any) => void;

export class KernelSelectionManager {
  private selectedAtomId: string | null = null;
  private listeners: Map<string, Set<Handler>> = new Map();

  constructor() {
    this.listeners.set('select', new Set());
    this.listeners.set('deselect', new Set());
    this.listeners.set('change', new Set());
  }

  getSelectedAtomId(): string | null {
    return this.selectedAtomId;
  }

  selectAtom(id: string): void {
    if (this.selectedAtomId === id) return;
    this.selectedAtomId = id;
    this.emit('select', id);
    this.emit('change', { selected: id });
  }

  deselect(): void {
    if (!this.selectedAtomId) return;
    const prev = this.selectedAtomId;
    this.selectedAtomId = null;
    this.emit('deselect', prev);
    this.emit('change', { selected: null });
  }

  on(event: 'select' | 'deselect' | 'change', cb: Handler): () => void {
    const set = this.listeners.get(event);
    set?.add(cb);
    return () => set?.delete(cb);
  }

  private emit(event: string, payload?: any): void {
    const set = this.listeners.get(event);
    if (!set) return;
    for (const h of Array.from(set)) {
      try {
        h(payload);
      } catch (e) {
        // swallow handler errors (kernels shouldn't crash tests)
        // eslint-disable-next-line no-console
        console.warn('KernelSelectionManager handler error', e);
      }
    }
  }

  reset(): void {
    this.selectedAtomId = null;
    this.listeners.forEach((s) => s.clear());
    this.listeners.set('select', new Set());
    this.listeners.set('deselect', new Set());
    this.listeners.set('change', new Set());
  }
}

// Export singleton instance
export const kernelSelectionManager = new KernelSelectionManager();

// Also export class for testing
export default KernelSelectionManager;

