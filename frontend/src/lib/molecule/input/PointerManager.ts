/**
 * PointerManager - Unified pointer event handling
 * 
 * Phase 3: Event & Input System Rewrite
 * 
 * Handles:
 * - Pointer down/move/up events
 * - Drag detection and gestures
 * - Bond creation gestures (drag from atom to atom)
 * - Click vs drag distinction
 * - Multi-touch support (future)
 */

export interface PointerEvent {
  type: 'down' | 'move' | 'up' | 'cancel'
  x: number
  y: number
  button: number
  buttons: number
  timestamp: number
  target?: EventTarget | null
}

export interface DragState {
  isDragging: boolean
  startX: number
  startY: number
  currentX: number
  currentY: number
  deltaX: number
  deltaY: number
  distance: number
  startTime: number
}

export interface PointerState {
  isDown: boolean
  isDragging: boolean
  dragState: DragState | null
  lastEvent: PointerEvent | null
}

// Minimum distance to consider a drag (in pixels)
const DRAG_THRESHOLD = 5

// Maximum time for a click (in milliseconds)
const CLICK_MAX_TIME = 300

export class PointerManager {
  private state: PointerState = {
    isDown: false,
    isDragging: false,
    dragState: null,
    lastEvent: null,
  }

  private listeners: Map<string, Set<(event: PointerEvent) => void>> = new Map()
  private dragListeners: Set<(state: DragState) => void> = new Set()

  /**
   * Process a pointer event
   */
  processEvent(nativeEvent: MouseEvent | TouchEvent, type: 'down' | 'move' | 'up' | 'cancel'): void {
    const point = this.getEventPoint(nativeEvent)
    if (!point) return

    const event: PointerEvent = {
      type,
      x: point.x,
      y: point.y,
      button: 'button' in nativeEvent ? nativeEvent.button : 0,
      buttons: 'buttons' in nativeEvent ? nativeEvent.buttons : 0,
      timestamp: Date.now(),
      target: nativeEvent.target,
    }

    this.handleEvent(event)
  }

  /**
   * Get point coordinates from event
   */
  private getEventPoint(event: MouseEvent | TouchEvent): { x: number; y: number } | null {
    if ('touches' in event) {
      // Touch event
      if (event.touches.length > 0) {
        return { x: event.touches[0].clientX, y: event.touches[0].clientY }
      }
      if (event.changedTouches.length > 0) {
        return { x: event.changedTouches[0].clientX, y: event.changedTouches[0].clientY }
      }
      return null
    } else {
      // Mouse event
      return { x: event.clientX, y: event.clientY }
    }
  }

  /**
   * Handle pointer event
   */
  private handleEvent(event: PointerEvent): void {
    this.state.lastEvent = event

    switch (event.type) {
      case 'down':
        this.handleDown(event)
        break
      case 'move':
        this.handleMove(event)
        break
      case 'up':
      case 'cancel':
        this.handleUp(event)
        break
    }

    // Notify listeners
    const listeners = this.listeners.get(event.type)
    if (listeners) {
      listeners.forEach(listener => listener(event))
    }
  }

  /**
   * Handle pointer down
   */
  private handleDown(event: PointerEvent): void {
    this.state.isDown = true
    this.state.isDragging = false
    this.state.dragState = {
      isDragging: false,
      startX: event.x,
      startY: event.y,
      currentX: event.x,
      currentY: event.y,
      deltaX: 0,
      deltaY: 0,
      distance: 0,
      startTime: event.timestamp,
    }
  }

  /**
   * Handle pointer move
   */
  private handleMove(event: PointerEvent): void {
    if (!this.state.isDown || !this.state.dragState) return

    const dragState = this.state.dragState
    dragState.currentX = event.x
    dragState.currentY = event.y
    dragState.deltaX = event.x - dragState.startX
    dragState.deltaY = event.y - dragState.startY
    dragState.distance = Math.sqrt(
      dragState.deltaX * dragState.deltaX + dragState.deltaY * dragState.deltaY
    )

    // Check if we've crossed the drag threshold
    if (!this.state.isDragging && dragState.distance > DRAG_THRESHOLD) {
      this.state.isDragging = true
      dragState.isDragging = true
    }

    if (this.state.isDragging) {
      // Notify drag listeners
      this.dragListeners.forEach(listener => listener(dragState))
    }
  }

  /**
   * Handle pointer up
   */
  private handleUp(event: PointerEvent): void {
    const wasDragging = this.state.isDragging
    const dragState = this.state.dragState

    this.state.isDown = false
    this.state.isDragging = false

    // Determine if this was a click (not a drag)
    const isClick = !wasDragging && dragState && 
      (event.timestamp - dragState.startTime) < CLICK_MAX_TIME

    if (dragState) {
      dragState.isDragging = false
    }

    this.state.dragState = null

    // Emit click event if it was a click
    if (isClick && dragState) {
      const clickEvent: PointerEvent = {
        ...event,
        type: 'down', // Treat as click
      }
      const clickListeners = this.listeners.get('click')
      if (clickListeners) {
        clickListeners.forEach(listener => listener(clickEvent))
      }
    }
  }

  /**
   * Add event listener
   */
  on(eventType: string, listener: (event: PointerEvent) => void): () => void {
    if (!this.listeners.has(eventType)) {
      this.listeners.set(eventType, new Set())
    }
    this.listeners.get(eventType)!.add(listener)

    // Return unsubscribe function
    return () => {
      const listeners = this.listeners.get(eventType)
      if (listeners) {
        listeners.delete(listener)
      }
    }
  }

  /**
   * Add drag listener
   */
  onDrag(listener: (state: DragState) => void): () => void {
    this.dragListeners.add(listener)
    return () => {
      this.dragListeners.delete(listener)
    }
  }

  /**
   * Get current state
   */
  getState(): PointerState {
    return { ...this.state }
  }

  /**
   * Check if currently dragging
   */
  isDragging(): boolean {
    return this.state.isDragging
  }

  /**
   * Get current drag state
   */
  getDragState(): DragState | null {
    return this.state.dragState ? { ...this.state.dragState } : null
  }

  /**
   * Cancel current operation
   */
  cancel(): void {
    if (this.state.isDown) {
      this.handleUp({
        type: 'cancel',
        x: this.state.lastEvent?.x ?? 0,
        y: this.state.lastEvent?.y ?? 0,
        button: 0,
        buttons: 0,
        timestamp: Date.now(),
      })
    }
  }

  /**
   * Reset state
   */
  reset(): void {
    this.state = {
      isDown: false,
      isDragging: false,
      dragState: null,
      lastEvent: null,
    }
  }
}

