/**
 * Bond class - represents a bond between two atoms
 * 
 * Immutable data structure with validation.
 */

import type { Bond } from './types'

export class BondImpl {
  private _data: Bond

  constructor(data: Bond) {
    this._data = { ...data }
    this.validate()
  }

  get id(): string {
    return this._data.id
  }

  get atom1(): string {
    return this._data.atom1
  }

  get atom2(): string {
    return this._data.atom2
  }

  get order(): number {
    return this._data.order
  }

  get type(): 'single' | 'double' | 'triple' | 'aromatic' {
    return this._data.type ?? this.inferType()
  }

  get stereo(): string | undefined {
    return this._data.stereo
  }

  get metadata(): Record<string, any> {
    return { ...this._data.metadata }
  }

  /**
   * Check if bond connects a specific atom
   */
  connects(atomId: string): boolean {
    return this._data.atom1 === atomId || this._data.atom2 === atomId
  }

  /**
   * Get the other atom in the bond
   */
  getOtherAtom(atomId: string): string | null {
    if (this._data.atom1 === atomId) return this._data.atom2
    if (this._data.atom2 === atomId) return this._data.atom1
    return null
  }

  /**
   * Create a new bond with updated properties
   */
  update(updates: Partial<Bond>): BondImpl {
    return new BondImpl({ ...this._data, ...updates })
  }

  /**
   * Convert to plain object
   */
  toJSON(): Bond {
    return { ...this._data }
  }

  /**
   * Validate bond data
   */
  private validate(): void {
    if (!this._data.id) {
      throw new Error('Bond must have an id')
    }
    if (!this._data.atom1 || !this._data.atom2) {
      throw new Error('Bond must connect two atoms')
    }
    if (this._data.atom1 === this._data.atom2) {
      throw new Error('Bond cannot connect atom to itself')
    }
    if (![1, 2, 3, 1.5].includes(this._data.order)) {
      throw new Error(`Invalid bond order: ${this._data.order}`)
    }
  }

  /**
   * Infer bond type from order
   */
  private inferType(): 'single' | 'double' | 'triple' | 'aromatic' {
    if (this._data.order === 1.5) return 'aromatic'
    if (this._data.order === 1) return 'single'
    if (this._data.order === 2) return 'double'
    if (this._data.order === 3) return 'triple'
    return 'single'
  }
}

