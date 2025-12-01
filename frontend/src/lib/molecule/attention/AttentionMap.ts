/**
 * Attention Map utilities
 * 
 * Phase 14: Attention Map Overlay (ML Explainability)
 * 
 * Handles attention data normalization and mapping to atoms/bonds.
 */

import type { Molecule } from '../Molecule'
import type { AttentionData } from '../prediction'

export interface AtomAttentionMap {
  [atomId: string]: number
}

export interface BondAttentionMap {
  [bondId: string]: number
}

export interface AttentionMap {
  atoms: AtomAttentionMap
  bonds: BondAttentionMap
  minValue: number
  maxValue: number
}

/**
 * Normalize attention values to [0, 1] range
 */
export function normalizeAttention(values: number[]): number[] {
  if (values.length === 0) return []
  
  const min = Math.min(...values)
  const max = Math.max(...values)
  const range = max - min
  
  if (range === 0) {
    return values.map(() => 0.5) // All same value -> middle
  }
  
  return values.map(v => (v - min) / range)
}

/**
 * Map attention data to molecule atoms and bonds
 */
export function mapAttentionToMolecule(
  molecule: Molecule,
  attentionData: AttentionData
): AttentionMap {
  const atoms = molecule.getAtoms()
  const bonds = molecule.getBonds()
  
  const atomMap: AtomAttentionMap = {}
  const bondMap: BondAttentionMap = {}
  
  // Map node importance to atoms
  if (attentionData.nodeImportance && attentionData.nodeImportance.length > 0) {
    const normalized = normalizeAttention(attentionData.nodeImportance)
    
    // Map by index (assuming order matches)
    atoms.forEach((atom, idx) => {
      if (idx < normalized.length) {
        atomMap[atom.id] = normalized[idx]
      }
    })
  }
  
  // Map edge attentions to bonds
  if (attentionData.edgeAttentions && attentionData.edgeIndex && attentionData.edgeIndex.length > 0) {
    const normalized = normalizeAttention(attentionData.edgeAttentions)
    const edgeIndex = attentionData.edgeIndex
    
    // Create atom ID to index mapping
    const atomIdToIndex = new Map<string, number>()
    atoms.forEach((atom, idx) => {
      atomIdToIndex.set(atom.id, idx)
    })
    
    // Map edges to bonds
    bonds.forEach((bond) => {
      const atom1Idx = atomIdToIndex.get(bond.atom1)
      const atom2Idx = atomIdToIndex.get(bond.atom2)
      
      if (atom1Idx !== undefined && atom2Idx !== undefined) {
        // Find matching edge in edgeIndex
        for (let e = 0; e < edgeIndex[0].length && e < normalized.length; e++) {
          const u = edgeIndex[0][e]
          const v = edgeIndex[1][e]
          
          if (
            (u === atom1Idx && v === atom2Idx) ||
            (u === atom2Idx && v === atom1Idx)
          ) {
            bondMap[bond.id] = normalized[e]
            break
          }
        }
      }
    })
  }
  
  // Calculate min/max for legend
  const allValues = [
    ...Object.values(atomMap),
    ...Object.values(bondMap),
  ]
  
  const minValue = allValues.length > 0 ? Math.min(...allValues) : 0
  const maxValue = allValues.length > 0 ? Math.max(...allValues) : 1
  
  return {
    atoms: atomMap,
    bonds: bondMap,
    minValue,
    maxValue,
  }
}

/**
 * Get color for attention value (0-1 range)
 * Uses a gradient from blue (low) to red (high)
 */
export function getAttentionColor(value: number): string {
  // Clamp to [0, 1]
  const clamped = Math.max(0, Math.min(1, value))
  
  // Blue (low) -> Cyan -> Yellow -> Red (high)
  if (clamped < 0.33) {
    // Blue to Cyan
    const t = clamped / 0.33
    const r = Math.floor(0)
    const g = Math.floor(t * 255)
    const b = Math.floor(255)
    return `rgb(${r}, ${g}, ${b})`
  } else if (clamped < 0.66) {
    // Cyan to Yellow
    const t = (clamped - 0.33) / 0.33
    const r = Math.floor(t * 255)
    const g = Math.floor(255)
    const b = Math.floor((1 - t) * 255)
    return `rgb(${r}, ${g}, ${b})`
  } else {
    // Yellow to Red
    const t = (clamped - 0.66) / 0.34
    const r = Math.floor(255)
    const g = Math.floor((1 - t) * 255)
    const b = Math.floor(0)
    return `rgb(${r}, ${g}, ${b})`
  }
}

/**
 * Get opacity for attention overlay (0-1 range)
 */
export function getAttentionOpacity(value: number, baseOpacity: number = 0.6): number {
  const clamped = Math.max(0, Math.min(1, value))
  return baseOpacity * (0.3 + 0.7 * clamped) // Fade from 30% to 100% of base opacity
}

