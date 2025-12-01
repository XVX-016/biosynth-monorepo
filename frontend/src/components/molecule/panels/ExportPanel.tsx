/**
 * ExportPanel - Export and save molecule
 * 
 * Phase 13: Save / Load / Export System
 * 
 * Features:
 * - Export as SMILES, MolBlock, JSON
 * - Export as SVG, PNG
 * - Copy to clipboard
 * - Save to LocalStorage
 */

import React, { useState, useRef } from 'react'
import { useEditorContext } from '../EditorContext'
import { toSMILES, toMolBlock, exportMoleculeAsJSON, exportMoleculeAsSVG } from '@/lib/molecule/export'
import { saveToLocalStorage, clearLocalStorage } from '@/lib/molecule/storage/autosave'
import { exportCanvasAsPNG } from '@/lib/molecule/export/image'

export function ExportPanel() {
  const { molecule } = useEditorContext()
  const [exporting, setExporting] = useState<string | null>(null)
  const fileInputRef = useRef<HTMLInputElement>(null)

  const isEmpty = molecule.isEmpty()

  // Copy SMILES to clipboard
  const handleCopySMILES = async () => {
    if (isEmpty) {
      alert('Molecule is empty')
      return
    }

    try {
      setExporting('smiles')
      const smiles = await toSMILES(molecule, true)
      await navigator.clipboard.writeText(smiles)
      alert('SMILES copied to clipboard!')
    } catch (error) {
      console.error('Error copying SMILES:', error)
      alert('Failed to copy SMILES')
    } finally {
      setExporting(null)
    }
  }

  // Export as SMILES file
  const handleExportSMILES = async () => {
    if (isEmpty) {
      alert('Molecule is empty')
      return
    }

    try {
      setExporting('smiles-file')
      const smiles = await toSMILES(molecule, true)
      const blob = new Blob([smiles], { type: 'text/plain' })
      const url = URL.createObjectURL(blob)
      const link = document.createElement('a')
      link.href = url
      link.download = 'molecule.smi'
      document.body.appendChild(link)
      link.click()
      document.body.removeChild(link)
      URL.revokeObjectURL(url)
    } catch (error) {
      console.error('Error exporting SMILES:', error)
      alert('Failed to export SMILES')
    } finally {
      setExporting(null)
    }
  }

  // Export as MolBlock
  const handleExportMolBlock = async () => {
    if (isEmpty) {
      alert('Molecule is empty')
      return
    }

    try {
      setExporting('molblock')
      const molblock = await toMolBlock(molecule)
      const blob = new Blob([molblock], { type: 'text/plain' })
      const url = URL.createObjectURL(blob)
      const link = document.createElement('a')
      link.href = url
      link.download = 'molecule.mol'
      document.body.appendChild(link)
      link.click()
      document.body.removeChild(link)
      URL.revokeObjectURL(url)
    } catch (error) {
      console.error('Error exporting MolBlock:', error)
      alert('Failed to export MolBlock')
    } finally {
      setExporting(null)
    }
  }

  // Export as JSON
  const handleExportJSON = () => {
    if (isEmpty) {
      alert('Molecule is empty')
      return
    }

    try {
      exportMoleculeAsJSON(molecule, 'molecule.json')
    } catch (error) {
      console.error('Error exporting JSON:', error)
      alert('Failed to export JSON')
    }
  }

  // Export as SVG
  const handleExportSVG = () => {
    if (isEmpty) {
      alert('Molecule is empty')
      return
    }

    try {
      exportMoleculeAsSVG(molecule, 'molecule.svg')
    } catch (error) {
      console.error('Error exporting SVG:', error)
      alert('Failed to export SVG')
    }
  }

  // Export as PNG (from canvas)
  const handleExportPNG = () => {
    if (isEmpty) {
      alert('Molecule is empty')
      return
    }

    // Get canvas from DOM (this is a workaround - ideally we'd pass canvas ref)
    const canvas = document.querySelector('canvas') as HTMLCanvasElement
    if (!canvas) {
      alert('Canvas not found. Please ensure the editor is visible.')
      return
    }

    try {
      exportCanvasAsPNG(canvas, 'molecule.png')
    } catch (error) {
      console.error('Error exporting PNG:', error)
      alert('Failed to export PNG')
    }
  }

  // Save to LocalStorage
  const handleSaveToLocalStorage = () => {
    if (isEmpty) {
      alert('Molecule is empty')
      return
    }

    try {
      saveToLocalStorage(molecule)
      alert('Molecule saved to LocalStorage!')
    } catch (error) {
      console.error('Error saving to LocalStorage:', error)
      alert('Failed to save to LocalStorage')
    }
  }

  // Clear LocalStorage
  const handleClearLocalStorage = () => {
    if (confirm('Clear autosave from LocalStorage?')) {
      clearLocalStorage()
      alert('LocalStorage cleared')
    }
  }

  return (
    <div className="bg-white border border-gray-200 rounded-lg p-4">
      <h3 className="text-sm font-semibold text-gray-700 mb-3">Export & Save</h3>

      <div className="space-y-2">
        {/* Text formats */}
        <div className="space-y-1">
          <label className="text-xs text-gray-500">Text Formats</label>
          <div className="grid grid-cols-2 gap-2">
            <button
              onClick={handleCopySMILES}
              disabled={isEmpty || exporting === 'smiles'}
              className="px-3 py-1.5 text-xs rounded border border-gray-300 bg-white text-gray-700 hover:bg-gray-50 disabled:opacity-50 disabled:cursor-not-allowed"
            >
              {exporting === 'smiles' ? 'Copying...' : 'Copy SMILES'}
            </button>
            <button
              onClick={handleExportSMILES}
              disabled={isEmpty || exporting === 'smiles-file'}
              className="px-3 py-1.5 text-xs rounded border border-gray-300 bg-white text-gray-700 hover:bg-gray-50 disabled:opacity-50 disabled:cursor-not-allowed"
            >
              {exporting === 'smiles-file' ? 'Exporting...' : 'Export SMILES'}
            </button>
            <button
              onClick={handleExportMolBlock}
              disabled={isEmpty || exporting === 'molblock'}
              className="px-3 py-1.5 text-xs rounded border border-gray-300 bg-white text-gray-700 hover:bg-gray-50 disabled:opacity-50 disabled:cursor-not-allowed"
            >
              {exporting === 'molblock' ? 'Exporting...' : 'Export MOL'}
            </button>
            <button
              onClick={handleExportJSON}
              disabled={isEmpty}
              className="px-3 py-1.5 text-xs rounded border border-gray-300 bg-white text-gray-700 hover:bg-gray-50 disabled:opacity-50 disabled:cursor-not-allowed"
            >
              Export JSON
            </button>
          </div>
        </div>

        {/* Image formats */}
        <div className="space-y-1">
          <label className="text-xs text-gray-500">Image Formats</label>
          <div className="grid grid-cols-2 gap-2">
            <button
              onClick={handleExportSVG}
              disabled={isEmpty}
              className="px-3 py-1.5 text-xs rounded border border-gray-300 bg-white text-gray-700 hover:bg-gray-50 disabled:opacity-50 disabled:cursor-not-allowed"
            >
              Export SVG
            </button>
            <button
              onClick={handleExportPNG}
              disabled={isEmpty}
              className="px-3 py-1.5 text-xs rounded border border-gray-300 bg-white text-gray-700 hover:bg-gray-50 disabled:opacity-50 disabled:cursor-not-allowed"
            >
              Export PNG
            </button>
          </div>
        </div>

        {/* LocalStorage */}
        <div className="space-y-1">
          <label className="text-xs text-gray-500">Local Storage</label>
          <div className="grid grid-cols-2 gap-2">
            <button
              onClick={handleSaveToLocalStorage}
              disabled={isEmpty}
              className="px-3 py-1.5 text-xs rounded border border-gray-300 bg-white text-gray-700 hover:bg-gray-50 disabled:opacity-50 disabled:cursor-not-allowed"
            >
              Save
            </button>
            <button
              onClick={handleClearLocalStorage}
              className="px-3 py-1.5 text-xs rounded border border-gray-300 bg-white text-gray-700 hover:bg-gray-50"
            >
              Clear
            </button>
          </div>
        </div>
      </div>
    </div>
  )
}

