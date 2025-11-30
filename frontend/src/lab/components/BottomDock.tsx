/**
 * BottomDock - Context-aware bottom panel
 * 
 * Shows different content based on selection:
 * - Precision editor when atom selected
 * - Bond inspector when bond selected
 * - Spectroscopy module when requested
 * - AI model output
 */

import React, { useState } from 'react'
import { toolController } from '../engines/ToolController'
import { moleculeEngine } from '../engines/MoleculeStateEngine'
import { ValidationEngine, type ValidationResult } from '../engines/ValidationEngine'

interface BottomDockProps {
  selectedAtomId?: string | null
  selectedBondId?: string | null
  validationResult?: ValidationResult | null
  mlPredictions?: Record<string, number> | null
}

export default function BottomDock({
  selectedAtomId,
  selectedBondId,
  validationResult,
  mlPredictions,
}: BottomDockProps) {
  const [activeTab, setActiveTab] = useState<'precision' | 'validation' | 'ml'>('precision')

  // Get selected atom data
  const selectedAtom = selectedAtomId ? moleculeEngine.getAtom(selectedAtomId) : null

  // Render precision editor
  const renderPrecisionEditor = () => {
    if (!selectedAtom) {
      return (
        <div className="p-4 text-center text-gray-500">
          Select an atom to edit its properties
        </div>
      )
    }

    return (
      <div className="p-4 space-y-4">
        <h3 className="text-lg font-semibold text-gray-800">Precision Editor</h3>
        
        <div className="grid grid-cols-2 gap-4">
          <div>
            <label className="block text-sm font-medium text-gray-700 mb-1">
              Element
            </label>
            <input
              type="text"
              value={selectedAtom.element}
              readOnly
              className="w-full p-2 border border-gray-300 rounded bg-gray-50"
            />
          </div>
          
          <div>
            <label className="block text-sm font-medium text-gray-700 mb-1">
              Charge
            </label>
            <input
              type="number"
              value={selectedAtom.charge}
              onChange={(e) => {
                const atom = moleculeEngine.getAtom(selectedAtomId!)
                if (atom) {
                  atom.charge = parseInt(e.target.value) || 0
                }
              }}
              className="w-full p-2 border border-gray-300 rounded"
            />
          </div>
          
          <div>
            <label className="block text-sm font-medium text-gray-700 mb-1">
              X Position
            </label>
            <input
              type="number"
              value={selectedAtom.x.toFixed(2)}
              onChange={(e) => {
                const atom = moleculeEngine.getAtom(selectedAtomId!)
                if (atom) {
                  moleculeEngine.setAtomPosition(selectedAtomId!, parseFloat(e.target.value) || 0, atom.y, atom.z)
                }
              }}
              className="w-full p-2 border border-gray-300 rounded"
            />
          </div>
          
          <div>
            <label className="block text-sm font-medium text-gray-700 mb-1">
              Y Position
            </label>
            <input
              type="number"
              value={selectedAtom.y.toFixed(2)}
              onChange={(e) => {
                const atom = moleculeEngine.getAtom(selectedAtomId!)
                if (atom) {
                  moleculeEngine.setAtomPosition(selectedAtomId!, atom.x, parseFloat(e.target.value) || 0, atom.z)
                }
              }}
              className="w-full p-2 border border-gray-300 rounded"
            />
          </div>
        </div>
      </div>
    )
  }

  // Render validation panel
  const renderValidation = () => {
    if (!validationResult) {
      return (
        <div className="p-4 text-center text-gray-500">
          Run validation to see results
        </div>
      )
    }

    return (
      <div className="p-4 space-y-3">
        <h3 className="text-lg font-semibold text-gray-800">Validation Results</h3>
        
        <div className={`p-3 rounded ${
          validationResult.valid ? 'bg-green-50 border border-green-200' : 'bg-red-50 border border-red-200'
        }`}>
          <p className={`font-medium ${
            validationResult.valid ? 'text-green-800' : 'text-red-800'
          }`}>
            {validationResult.valid ? '✓ Structure is valid' : '✗ Structure has errors'}
          </p>
        </div>

        {validationResult.issues.length > 0 && (
          <div className="space-y-2">
            {validationResult.issues.map((issue, idx) => (
              <div
                key={idx}
                className={`p-2 rounded text-sm ${
                  issue.type === 'error' ? 'bg-red-100 text-red-800' :
                  issue.type === 'warning' ? 'bg-yellow-100 text-yellow-800' :
                  'bg-blue-100 text-blue-800'
                }`}
              >
                {issue.message}
              </div>
            ))}
          </div>
        )}
      </div>
    )
  }

  // Render ML predictions
  const renderMLPredictions = () => {
    if (!mlPredictions) {
      return (
        <div className="p-4 text-center text-gray-500">
          Run property prediction to see results
        </div>
      )
    }

    return (
      <div className="p-4 space-y-3">
        <h3 className="text-lg font-semibold text-gray-800">ML Predictions</h3>
        
        <div className="grid grid-cols-2 gap-4">
          {Object.entries(mlPredictions).map(([key, value]) => (
            <div key={key} className="p-3 bg-gray-50 rounded border border-gray-200">
              <p className="text-sm text-gray-600 capitalize">{key}</p>
              <p className="text-xl font-semibold text-gray-800">{value.toFixed(3)}</p>
            </div>
          ))}
        </div>
      </div>
    )
  }

  return (
    <div className="w-full h-64 bg-white border-t border-gray-300 flex flex-col">
      {/* Tabs */}
      <div className="flex border-b border-gray-300">
        <button
          onClick={() => setActiveTab('precision')}
          className={`px-4 py-2 text-sm font-medium ${
            activeTab === 'precision'
              ? 'text-blue-700 border-b-2 border-blue-700'
              : 'text-gray-600 hover:text-gray-800'
          }`}
        >
          Precision Editor
        </button>
        <button
          onClick={() => setActiveTab('validation')}
          className={`px-4 py-2 text-sm font-medium ${
            activeTab === 'validation'
              ? 'text-blue-700 border-b-2 border-blue-700'
              : 'text-gray-600 hover:text-gray-800'
          }`}
        >
          Validation
        </button>
        <button
          onClick={() => setActiveTab('ml')}
          className={`px-4 py-2 text-sm font-medium ${
            activeTab === 'ml'
              ? 'text-blue-700 border-b-2 border-blue-700'
              : 'text-gray-600 hover:text-gray-800'
          }`}
        >
          ML Predictions
        </button>
      </div>

      {/* Content */}
      <div className="flex-1 overflow-y-auto">
        {activeTab === 'precision' && renderPrecisionEditor()}
        {activeTab === 'validation' && renderValidation()}
        {activeTab === 'ml' && renderMLPredictions()}
      </div>
    </div>
  )
}

