/**
 * InspectorPanel - Atom/Bond property inspector
 * 
 * Phase 10: Lab Page UI Rebuild
 * Phase 12: Will be expanded with editable fields
 */

import React from 'react'
import { useEditorContext } from '../EditorContext'

export function InspectorPanel() {
  const { molecule, selectedAtomId, selectedBondId } = useEditorContext()

  const selectedAtom = selectedAtomId ? molecule.getAtom(selectedAtomId) : null
  const selectedBond = selectedBondId ? molecule.getBond(selectedBondId) : null

  if (!selectedAtom && !selectedBond) {
    return (
      <div className="bg-white border border-gray-200 rounded-lg p-4">
        <h3 className="text-sm font-semibold text-gray-700 mb-2">Inspector</h3>
        <p className="text-xs text-gray-500">Select an atom or bond to inspect</p>
      </div>
    )
  }

  return (
    <div className="bg-white border border-gray-200 rounded-lg p-4">
      <h3 className="text-sm font-semibold text-gray-700 mb-3">Inspector</h3>
      
      {selectedAtom && (
        <div className="space-y-2">
          <div>
            <label className="text-xs text-gray-500">Element</label>
            <div className="text-sm font-mono">{selectedAtom.element}</div>
          </div>
          <div>
            <label className="text-xs text-gray-500">Position</label>
            <div className="text-xs font-mono text-gray-600">
              ({selectedAtom.position[0].toFixed(2)}, {selectedAtom.position[1].toFixed(2)}, {selectedAtom.position[2].toFixed(2)})
            </div>
          </div>
          {selectedAtom.charge !== undefined && (
            <div>
              <label className="text-xs text-gray-500">Charge</label>
              <div className="text-sm">{selectedAtom.charge}</div>
            </div>
          )}
        </div>
      )}

      {selectedBond && (
        <div className="space-y-2 mt-4">
          <div>
            <label className="text-xs text-gray-500">Bond Order</label>
            <div className="text-sm font-mono">{selectedBond.order}</div>
          </div>
          <div>
            <label className="text-xs text-gray-500">Atoms</label>
            <div className="text-xs text-gray-600">
              {selectedBond.atom1} - {selectedBond.atom2}
            </div>
          </div>
        </div>
      )}
    </div>
  )
}

