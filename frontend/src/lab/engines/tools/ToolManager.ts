/**
 * ToolManager - Manages active tool and tool registration
 */

import { ITool } from './Tool.types';

export class ToolManager {
  private tools: Record<string, ITool> = {};
  private active: ITool | null = null;
  private listeners: ((tool: ITool | null) => void)[] = [];

  /**
   * Register a tool
   */
  register(tool: ITool): void {
    this.tools[tool.id] = tool;
  }

  /**
   * Set the active tool
   */
  setActive(id: string): void {
    if (this.active?.onDeselect) {
      this.active.onDeselect();
    }
    
    this.active = this.tools[id] || null;
    
    if (this.active?.onSelect) {
      this.active.onSelect();
    }
    
    this.notify(this.active);
  }

  /**
   * Get the active tool
   */
  getActive(): ITool | null {
    return this.active;
  }

  /**
   * Get tool by ID
   */
  getTool(id: string): ITool | undefined {
    return this.tools[id];
  }

  /**
   * Get all registered tools
   */
  getAllTools(): ITool[] {
    return Object.values(this.tools);
  }

  /**
   * Subscribe to tool changes
   */
  onChange(cb: (tool: ITool | null) => void): () => void {
    this.listeners.push(cb);
    return () => {
      const index = this.listeners.indexOf(cb);
      if (index > -1) {
        this.listeners.splice(index, 1);
      }
    };
  }

  /**
   * Notify listeners of tool change
   */
  private notify(tool: ITool | null): void {
    this.listeners.forEach((cb) => cb(tool));
  }
}

