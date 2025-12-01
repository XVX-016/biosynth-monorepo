/**
 * KeyboardManager - Unified keyboard event handling
 * 
 * Phase 3: Event & Input System Rewrite
 * 
 * Handles:
 * - Keyboard shortcuts
 * - Escape to cancel operations
 * - Delete/Backspace to remove selected items
 * - Arrow keys for navigation (future)
 * - Undo/Redo shortcuts (future)
 */

export interface KeyboardEvent {
  key: string
  code: string
  ctrlKey: boolean
  shiftKey: boolean
  altKey: boolean
  metaKey: boolean
  timestamp: number
  preventDefault?: () => void
  stopPropagation?: () => void
}

export type KeyboardShortcut = {
  key: string
  ctrl?: boolean
  shift?: boolean
  alt?: boolean
  meta?: boolean
  handler: (event: KeyboardEvent) => void
  description?: string
}

export class KeyboardManager {
  private shortcuts: Map<string, KeyboardShortcut> = new Map()
  private listeners: Set<(event: KeyboardEvent) => void> = new Set()
  private enabled: boolean = true

  constructor() {
    // Register default shortcuts
    this.registerDefaults()
  }

  /**
   * Register default keyboard shortcuts
   */
  private registerDefaults(): void {
    // Escape - cancel current operation
    this.register({
      key: 'Escape',
      handler: (event) => {
        this.emit('cancel', event)
      },
      description: 'Cancel current operation',
    })

    // Delete - delete selected item
    this.register({
      key: 'Delete',
      handler: (event) => {
        this.emit('delete', event)
      },
      description: 'Delete selected item',
    })

    // Backspace - delete selected item
    this.register({
      key: 'Backspace',
      handler: (event) => {
        this.emit('delete', event)
      },
      description: 'Delete selected item',
    })
  }

  /**
   * Register a keyboard shortcut
   */
  register(shortcut: KeyboardShortcut): () => void {
    const id = this.getShortcutId(shortcut)
    this.shortcuts.set(id, shortcut)

    // Return unsubscribe function
    return () => {
      this.shortcuts.delete(id)
    }
  }

  /**
   * Process keyboard event
   */
  processEvent(nativeEvent: KeyboardEvent): boolean {
    if (!this.enabled) return false

    const event: KeyboardEvent = {
      key: nativeEvent.key,
      code: nativeEvent.code,
      ctrlKey: nativeEvent.ctrlKey,
      shiftKey: nativeEvent.shiftKey,
      altKey: nativeEvent.altKey,
      metaKey: nativeEvent.metaKey,
      timestamp: Date.now(),
      preventDefault: nativeEvent.preventDefault,
      stopPropagation: nativeEvent.stopPropagation,
    }

    // Check for matching shortcut
    const shortcut = this.findMatchingShortcut(event)
    if (shortcut) {
      shortcut.handler(event)
      if (event.preventDefault) {
        event.preventDefault()
      }
      return true
    }

    // Notify general listeners
    this.listeners.forEach(listener => listener(event))
    return false
  }

  /**
   * Find matching shortcut for event
   */
  private findMatchingShortcut(event: KeyboardEvent): KeyboardShortcut | null {
    for (const shortcut of this.shortcuts.values()) {
      if (this.matchesShortcut(event, shortcut)) {
        return shortcut
      }
    }
    return null
  }

  /**
   * Check if event matches shortcut
   */
  private matchesShortcut(event: KeyboardEvent, shortcut: KeyboardShortcut): boolean {
    if (event.key !== shortcut.key) return false
    if (shortcut.ctrl !== undefined && event.ctrlKey !== shortcut.ctrl) return false
    if (shortcut.shift !== undefined && event.shiftKey !== shortcut.shift) return false
    if (shortcut.alt !== undefined && event.altKey !== shortcut.alt) return false
    if (shortcut.meta !== undefined && event.metaKey !== shortcut.meta) return false
    return true
  }

  /**
   * Get shortcut ID
   */
  private getShortcutId(shortcut: KeyboardShortcut): string {
    const parts = [
      shortcut.ctrl ? 'ctrl' : '',
      shortcut.shift ? 'shift' : '',
      shortcut.alt ? 'alt' : '',
      shortcut.meta ? 'meta' : '',
      shortcut.key.toLowerCase(),
    ].filter(Boolean)
    return parts.join('+')
  }

  /**
   * Add general keyboard listener
   */
  on(listener: (event: KeyboardEvent) => void): () => void {
    this.listeners.add(listener)
    return () => {
      this.listeners.delete(listener)
    }
  }

  /**
   * Emit custom event
   */
  private emit(eventType: string, event: KeyboardEvent): void {
    const listeners = this.listeners
    // Could add event-specific listeners here
    listeners.forEach(listener => listener(event))
  }

  /**
   * Enable/disable keyboard handling
   */
  setEnabled(enabled: boolean): void {
    this.enabled = enabled
  }

  /**
   * Get all registered shortcuts
   */
  getShortcuts(): KeyboardShortcut[] {
    return Array.from(this.shortcuts.values())
  }

  /**
   * Get shortcut description
   */
  getShortcutDescription(key: string): string | undefined {
    for (const shortcut of this.shortcuts.values()) {
      if (shortcut.key === key) {
        return shortcut.description
      }
    }
    return undefined
  }
}

