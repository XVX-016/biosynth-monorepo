import React from 'react'
import { useLabStore } from '../../store/labStore'

export default function Inspector() {
  const mol = useLabStore(s => s.molecule)
  const auto = useLabStore(s => s.autoBond)
  const setAuto = useLabStore(s => s.setAutoBond)
  const undo = useLabStore(s => s.undo)
  const redo = useLabStore(s => s.redo)
  const resetMolecule = useLabStore(s => s.resetMolecule)
  
  return (
    <div className="bg-white rounded-lg border border-gray-200 p-3 shadow-sm">
      <h4 className="text-xs font-semibold text-gray-700 mb-3 text-center">Inspector</h4>
      <div className="space-y-2 text-xs">
        <div className="flex justify-between">
          <span className="text-gray-600">Atoms:</span>
          <span className="font-semibold text-gray-900">{mol.atoms.length}</span>
        </div>
        <div className="flex justify-between">
          <span className="text-gray-600">Bonds:</span>
          <span className="font-semibold text-gray-900">{mol.bonds.length}</span>
        </div>
      </div>
      <div className="mt-3 pt-3 border-t border-gray-200">
        <label className="flex items-center gap-2 text-xs cursor-pointer">
          <input
            type="checkbox"
            checked={auto}
            onChange={(e) => setAuto(e.target.checked)}
            className="rounded"
          />
          <span className="text-gray-700">Auto-bond</span>
        </label>
      </div>
      <div className="mt-3 pt-3 border-t border-gray-200 grid grid-cols-3 gap-1">
        <button
          onClick={undo}
          className="px-2 py-1 text-[10px] rounded border border-gray-300 bg-white hover:bg-gray-50 transition-colors"
          title="Undo"
        >
          ↶
        </button>
        <button
          onClick={redo}
          className="px-2 py-1 text-[10px] rounded border border-gray-300 bg-white hover:bg-gray-50 transition-colors"
          title="Redo"
        >
          ↷
        </button>
        <button
          onClick={resetMolecule}
          className="px-2 py-1 text-[10px] rounded border border-red-300 bg-white hover:bg-red-50 text-red-600 transition-colors"
          title="Reset"
        >
          ×
        </button>
      </div>
    </div>
  )
}

