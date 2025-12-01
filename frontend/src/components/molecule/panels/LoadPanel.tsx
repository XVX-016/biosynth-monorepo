/**
 * LoadPanel - Load molecule from file or input
 * 
 * Phase 13: Save / Load / Export System
 * 
 * Features:
 * - Load from SMILES input
 * - Load from file (JSON, MOL, SMILES)
 * - Drag-and-drop support
 */

import React, { useState, useRef } from 'react'
import { useEditorContext } from '../EditorContext'
import { loadFromSMILES, loadFromFile, loadFromDropEvent } from '@/lib/molecule/load'
import { loadFromLocalStorage } from '@/lib/molecule/storage/autosave'

export function LoadPanel() {
  const { setMolecule } = useEditorContext()
  const [smilesInput, setSmilesInput] = useState('')
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const fileInputRef = useRef<HTMLInputElement>(null)

  // Load from SMILES input
  const handleLoadFromSMILES = async () => {
    if (!smilesInput.trim()) {
      setError('Please enter a SMILES string')
      return
    }

    try {
      setLoading(true)
      setError(null)
      const molecule = await loadFromSMILES(smilesInput.trim())
      setMolecule(molecule)
      setSmilesInput('')
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to load SMILES')
    } finally {
      setLoading(false)
    }
  }

  // Load from file
  const handleFileSelect = async (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0]
    if (!file) return

    try {
      setLoading(true)
      setError(null)
      const molecule = await loadFromFile(file)
      setMolecule(molecule)
      if (fileInputRef.current) {
        fileInputRef.current.value = ''
      }
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to load file')
    } finally {
      setLoading(false)
    }
  }

  // Handle drag and drop
  const handleDragOver = (event: React.DragEvent) => {
    event.preventDefault()
    event.stopPropagation()
  }

  const handleDrop = async (event: React.DragEvent) => {
    event.preventDefault()
    event.stopPropagation()

    try {
      setLoading(true)
      setError(null)
      const molecule = await loadFromDropEvent(event.nativeEvent)
      if (molecule) {
        setMolecule(molecule)
      }
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to load file')
    } finally {
      setLoading(false)
    }
  }

  // Load from LocalStorage
  const handleLoadFromLocalStorage = () => {
    try {
      const molecule = loadFromLocalStorage()
      if (molecule) {
        setMolecule(molecule)
        setError(null)
      } else {
        setError('No saved molecule found in LocalStorage')
      }
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to load from LocalStorage')
    }
  }

  return (
    <div className="bg-white border border-gray-200 rounded-lg p-4">
      <h3 className="text-sm font-semibold text-gray-700 mb-3">Load Molecule</h3>

      <div className="space-y-3">
        {/* SMILES input */}
        <div>
          <label className="text-xs text-gray-500 mb-1 block">SMILES String</label>
          <div className="flex gap-2">
            <input
              type="text"
              value={smilesInput}
              onChange={(e) => setSmilesInput(e.target.value)}
              onKeyPress={(e) => e.key === 'Enter' && handleLoadFromSMILES()}
              placeholder="e.g., CCO (ethanol)"
              className="flex-1 px-2 py-1.5 text-xs border border-gray-300 rounded focus:outline-none focus:ring-2 focus:ring-blue-500"
              disabled={loading}
            />
            <button
              onClick={handleLoadFromSMILES}
              disabled={loading || !smilesInput.trim()}
              className="px-3 py-1.5 text-xs rounded border border-gray-300 bg-white text-gray-700 hover:bg-gray-50 disabled:opacity-50 disabled:cursor-not-allowed"
            >
              Load
            </button>
          </div>
        </div>

        {/* File upload */}
        <div>
          <label className="text-xs text-gray-500 mb-1 block">File Upload</label>
          <input
            ref={fileInputRef}
            type="file"
            accept=".json,.mol,.sdf,.smi,.smiles"
            onChange={handleFileSelect}
            className="hidden"
            id="molecule-file-input"
          />
          <label
            htmlFor="molecule-file-input"
            className="block w-full px-3 py-1.5 text-xs rounded border border-gray-300 bg-white text-gray-700 hover:bg-gray-50 cursor-pointer text-center"
          >
            Choose File (JSON, MOL, SMILES)
          </label>
        </div>

        {/* Drag and drop */}
        <div
          onDragOver={handleDragOver}
          onDrop={handleDrop}
          className="border-2 border-dashed border-gray-300 rounded p-4 text-center hover:border-gray-400 transition-colors"
        >
          <p className="text-xs text-gray-500">
            Drag and drop a molecule file here
          </p>
          <p className="text-xs text-gray-400 mt-1">
            Supports: JSON, MOL, SDF, SMILES
          </p>
        </div>

        {/* Load from LocalStorage */}
        <div>
          <button
            onClick={handleLoadFromLocalStorage}
            className="w-full px-3 py-1.5 text-xs rounded border border-gray-300 bg-white text-gray-700 hover:bg-gray-50"
          >
            Load from Autosave
          </button>
        </div>

        {/* Error message */}
        {error && (
          <div className="text-xs text-red-600 bg-red-50 p-2 rounded">
            {error}
          </div>
        )}

        {/* Loading indicator */}
        {loading && (
          <div className="text-xs text-gray-500 text-center">
            Loading...
          </div>
        )}
      </div>
    </div>
  )
}

