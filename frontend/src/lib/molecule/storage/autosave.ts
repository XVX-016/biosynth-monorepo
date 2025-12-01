/**
 * Autosave functionality
 * 
 * Phase 13: Save / Load / Export System
 * 
 * Handles LocalStorage autosave and session restoration.
 */

import type { Molecule } from '../Molecule'
import { toJSON, fromJSON } from '../export/json'

const AUTOSAVE_KEY = 'molforge_autosave'
const AUTOSAVE_INTERVAL = 5000 // 5 seconds

/**
 * Save molecule to LocalStorage
 */
export function saveToLocalStorage(molecule: Molecule): void {
  try {
    const json = toJSON(molecule)
    localStorage.setItem(AUTOSAVE_KEY, json)
  } catch (error) {
    console.error('Failed to save to LocalStorage:', error)
  }
}

/**
 * Load molecule from LocalStorage
 */
export function loadFromLocalStorage(): Molecule | null {
  try {
    const json = localStorage.getItem(AUTOSAVE_KEY)
    if (!json) {
      return null
    }
    return fromJSON(json)
  } catch (error) {
    console.error('Failed to load from LocalStorage:', error)
    return null
  }
}

/**
 * Clear autosave from LocalStorage
 */
export function clearLocalStorage(): void {
  try {
    localStorage.removeItem(AUTOSAVE_KEY)
  } catch (error) {
    console.error('Failed to clear LocalStorage:', error)
  }
}

/**
 * Setup autosave interval
 */
export function setupAutosave(
  molecule: Molecule,
  interval: number = AUTOSAVE_INTERVAL
): () => void {
  // Save immediately
  saveToLocalStorage(molecule)
  
  // Setup interval
  const timer = setInterval(() => {
    saveToLocalStorage(molecule)
  }, interval)
  
  // Return cleanup function
  return () => {
    clearInterval(timer)
  }
}

