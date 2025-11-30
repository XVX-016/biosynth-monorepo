/**
 * ToolController - Manages active tool mode and tool-specific event handling
 */

export enum ToolMode {
  cursor = "cursor",
  addAtom = "addAtom",
  addBond = "addBond",
  erase = "erase",
  select = "select",
  measure = "measure"
}

export interface Point {
  x: number
  y: number
}

export class ToolController {
  activeTool: ToolMode = ToolMode.cursor
  selectedAtom: string | null = null
  bondPreviewStart: string | null = null // For bond drawing preview

  /**
   * Set the active tool
   */
  setTool(tool: ToolMode): void {
    this.activeTool = tool
    // Reset tool-specific state
    this.selectedAtom = null
    this.bondPreviewStart = null
  }

  /**
   * Handle mouse down event
   * Routes to appropriate tool handler
   */
  handleMouseDown(pos: Point, hitAtomId: string | null): void {
    switch (this.activeTool) {
      case ToolMode.cursor:
        this.selectedAtom = hitAtomId
        break
      
      case ToolMode.addAtom:
        // Handled by parent component (needs element selection)
        break
      
      case ToolMode.addBond:
        if (hitAtomId) {
          if (this.bondPreviewStart === null) {
            // First click - start bond
            this.bondPreviewStart = hitAtomId
          } else if (this.bondPreviewStart !== hitAtomId) {
            // Second click - complete bond
            // This will be handled by parent component
            return
          }
        }
        break
      
      case ToolMode.erase:
        if (hitAtomId) {
          // Will be handled by parent component
          return
        }
        break
      
      case ToolMode.select:
        this.selectedAtom = hitAtomId
        break
    }
  }

  /**
   * Handle mouse move event
   * Used for previews (bond drawing, etc.)
   */
  handleMouseMove(pos: Point, hitAtomId: string | null): void {
    // Preview logic can be added here
    // For example, bond preview line
  }

  /**
   * Handle mouse up event
   * Finalizes actions
   */
  handleMouseUp(pos: Point, hitAtomId: string | null): void {
    // Most actions are finalized on mouse down
    // This can be used for drag operations
  }

  /**
   * Reset tool state
   */
  reset(): void {
    this.selectedAtom = null
    this.bondPreviewStart = null
  }
}

// Singleton instance
export const toolController = new ToolController()

