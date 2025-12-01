/**
 * EditorPanel - Main molecule editor panel
 * 
 * Phase 10: Lab Page UI Rebuild
 */

import React, { useState, useEffect } from 'react'
import { MoleculeEditor } from '../MoleculeEditor'
import { ThreeDViewer } from '../ThreeDViewer'
import { useEditorContext } from '../EditorContext'
import { generate3DCoordinates } from '@/lib/molecule/3d'

export function EditorPanel() {
  const { molecule, setMolecule, selectedAtomId, setSelectedAtomId, tool } = useEditorContext()
  const [viewMode, setViewMode] = useState<'2d' | '3d'>('2d')
  const [threeDMolecule, setThreeDMolecule] = useState<typeof molecule | null>(null)
  const [loading3D, setLoading3D] = useState(false)

  // Generate 3D coordinates when switching to 3D view
  const handleViewModeChange = async (mode: '2d' | '3d') => {
    setViewMode(mode)
    
    if (mode === '3d' && !threeDMolecule && !molecule.isEmpty()) {
      setLoading3D(true)
      try {
        const mol3d = await generate3DCoordinates(molecule)
        setThreeDMolecule(mol3d)
      } catch (error) {
        console.error('Failed to generate 3D coordinates:', error)
      } finally {
        setLoading3D(false)
      }
    }
  }

  return (
    <div className="flex flex-col h-full bg-white border border-gray-200 rounded-lg overflow-hidden">
      {/* View mode tabs */}
      <div className="flex border-b border-gray-200">
        <button
          onClick={() => handleViewModeChange('2d')}
          className={`flex-1 px-4 py-2 text-sm font-medium ${
            viewMode === '2d'
              ? 'bg-blue-50 text-blue-700 border-b-2 border-blue-700'
              : 'text-gray-600 hover:bg-gray-50'
          }`}
        >
          2D Editor
        </button>
        <button
          onClick={() => handleViewModeChange('3d')}
          className={`flex-1 px-4 py-2 text-sm font-medium ${
            viewMode === '3d'
              ? 'bg-blue-50 text-blue-700 border-b-2 border-blue-700'
              : 'text-gray-600 hover:bg-gray-50'
          }`}
        >
          3D Viewer
        </button>
      </div>

      {/* Editor content */}
      <div className="flex-1 relative">
        {viewMode === '2d' ? (
          <MoleculeEditor
            width={800}
            height={600}
            initialMolecule={molecule}
            tool={tool}
            onMoleculeChange={setMolecule}
            onAtomSelect={setSelectedAtomId}
          />
        ) : (
          <div className="w-full h-full">
            {loading3D ? (
              <div className="flex items-center justify-center h-full">
                <p className="text-sm text-gray-500">Generating 3D coordinates...</p>
              </div>
            ) : (
              <ThreeDViewer
                molecule={threeDMolecule || molecule}
                selectedAtomId={selectedAtomId}
                width={800}
                height={600}
              />
            )}
          </div>
        )}
      </div>
    </div>
  )
}

