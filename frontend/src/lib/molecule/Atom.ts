/**
 * Atom class - represents a single atom in a molecule
 * 
 * Immutable data structure with validation.
 */

import type { Atom, ElementInfo } from './types'
import { ELEMENT_DATA } from './constants'

export class AtomImpl {
  private _data: Atom

  constructor(data: Atom) {
    this._data = { ...data }
    this.validate()
  }

  get id(): string {
    return this._data.id
  }

  get element(): string {
    return this._data.element
  }

  get position(): [number, number, number] {
    return [...this._data.position] as [number, number, number]
  }

  get charge(): number {
    return this._data.charge ?? 0
  }

  get formalCharge(): number {
    return this._data.formalCharge ?? 0
  }

  get hybridization(): string | undefined {
    return this._data.hybridization
  }

  get aromatic(): boolean {
    return this._data.aromatic ?? false
  }

  get inRing(): boolean {
    return this._data.inRing ?? false
  }

  get valence(): number {
    return this._data.valence ?? 0
  }

  get maxValence(): number {
    return this._data.maxValence ?? this.getElementInfo().maxValence
  }

  get metadata(): Record<string, any> {
    return { ...this._data.metadata }
  }

  /**
   * Get element information
   */
  getElementInfo(): ElementInfo {
    const info = ELEMENT_DATA[this._data.element]
    if (!info) {
      throw new Error(`Unknown element: ${this._data.element}`)
    }
    return info
  }

  /**
   * Check if atom can accept more bonds
   */
  canAcceptBond(order: number = 1): boolean {
    return this.valence + order <= this.maxValence
  }

  /**
   * Create a new atom with updated properties
   */
  update(updates: Partial<Atom>): AtomImpl {
    return new AtomImpl({ ...this._data, ...updates })
  }

  /**
   * Convert to plain object
   */
  toJSON(): Atom {
    return { ...this._data }
  }

  /**
   * Validate atom data
   */
  private validate(): void {
    if (!this._data.id) {
      throw new Error('Atom must have an id')
    }
    if (!this._data.element) {
      throw new Error('Atom must have an element')
    }
    if (!ELEMENT_DATA[this._data.element]) {
      throw new Error(`Invalid element: ${this._data.element}`)
    }
    if (!this._data.position || this._data.position.length !== 3) {
      throw new Error('Atom must have valid 3D position')
    }
  }

  /**
   * Calculate distance to another atom
   */
  distanceTo(other: AtomImpl): number {
    const [x1, y1, z1] = this.position
    const [x2, y2, z2] = other.position
    return Math.sqrt(
      Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2) + Math.pow(z2 - z1, 2)
    )
  }
}

