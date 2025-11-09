import { useMoleculeStore } from '../../store/moleculeStore'

/**
 * SelectionManager - Tracks atom selection, hover, and drag state
 * Simple event emitter pattern for 3D interactions
 */
class SelectionManager {
  private hoveredAtomId: string | null = null
  private selectedAtomId: string | null = null
  private draggingAtomId: string | null = null
  private listeners: Map<string, Set<() => void>> = new Map()

  /**
   * Subscribe to state changes
   */
  on(event: 'hover' | 'select' | 'drag', callback: () => void): () => void {
    if (!this.listeners.has(event)) {
      this.listeners.set(event, new Set())
    }
    this.listeners.get(event)!.add(callback)

    // Return unsubscribe function
    return () => {
      this.listeners.get(event)?.delete(callback)
    }
  }

  /**
   * Emit event to all listeners
   */
  private emit(event: 'hover' | 'select' | 'drag'): void {
    this.listeners.get(event)?.forEach((callback) => callback())
  }

  /**
   * Handle atom hover
   */
  onHover(id: string | null): void {
    if (this.hoveredAtomId !== id) {
      this.hoveredAtomId = id
      this.emit('hover')
    }
  }

  /**
   * Handle atom selection
   */
  onSelect(id: string | null): void {
    if (this.selectedAtomId !== id) {
      this.selectedAtomId = id
      // Update store
      useMoleculeStore.getState().selectAtom(id)
      this.emit('select')
    }
  }

  /**
   * Start dragging an atom
   */
  startDrag(id: string): void {
    if (this.draggingAtomId !== id) {
      this.draggingAtomId = id
      this.emit('drag')
    }
  }

  /**
   * End dragging
   */
  endDrag(): void {
    if (this.draggingAtomId !== null) {
      this.draggingAtomId = null
      this.emit('drag')
    }
  }

  /**
   * Get current hovered atom ID
   */
  getHoveredAtomId(): string | null {
    return this.hoveredAtomId
  }

  /**
   * Get current selected atom ID
   */
  getSelectedAtomId(): string | null {
    return this.selectedAtomId
  }

  /**
   * Get current dragging atom ID
   */
  getDraggingAtomId(): string | null {
    return this.draggingAtomId
  }

  /**
   * Check if atom is hovered
   */
  isHovered(id: string): boolean {
    return this.hoveredAtomId === id
  }

  /**
   * Check if atom is selected
   */
  isSelected(id: string): boolean {
    return this.selectedAtomId === id
  }

  /**
   * Check if atom is being dragged
   */
  isDragging(id: string): boolean {
    return this.draggingAtomId === id
  }

  /**
   * Reset all state
   */
  reset(): void {
    this.hoveredAtomId = null
    this.selectedAtomId = null
    this.draggingAtomId = null
    this.emit('hover')
    this.emit('select')
    this.emit('drag')
  }
}

// Singleton instance
export const selectionManager = new SelectionManager()

