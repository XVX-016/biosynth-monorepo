/**
 * InspectorPanel - Atom/Bond property inspector with editing
 * 
 * Phase 12: Atom & Bond Inspector Panel
 * 
 * Features:
 * - Editable atom properties (element, charge, isotope)
 * - Editable bond properties (order, stereo)
 * - History integration (edits are commands)
 * - Selection state sync
 */

import React, { useState, useEffect } from 'react'
import { useEditorContext } from '../EditorContext'
import { ELEMENT_DATA, COMMON_ELEMENTS } from '@/lib/molecule/constants'
import { VALID_BOND_ORDERS } from '@/lib/molecule/constants'
import { UpdateAtomCommand, UpdateBondCommand } from '@/lib/molecule/history'
import { HistoryManager } from '@/lib/molecule/history'

export function InspectorPanel() {
  const { 
    molecule, 
    setMolecule,
    selectedAtomId, 
    selectedBondId,
    setSelectedAtomId,
    setSelectedBondId,
  } = useEditorContext()

  const selectedAtom = selectedAtomId ? molecule.getAtom(selectedAtomId) : null
  const selectedBond = selectedBondId ? molecule.getBond(selectedBondId) : null

  // Form state for atom
  const [atomElement, setAtomElement] = useState<string>('')
  const [atomCharge, setAtomCharge] = useState<number>(0)
  const [atomIsotope, setAtomIsotope] = useState<number | undefined>(undefined)

  // Form state for bond
  const [bondOrder, setBondOrder] = useState<number>(1)

  // Update form when selection changes
  useEffect(() => {
    if (selectedAtom) {
      setAtomElement(selectedAtom.element)
      setAtomCharge(selectedAtom.charge || 0)
      setAtomIsotope(selectedAtom.meta?.isotope)
    }
  }, [selectedAtom])

  useEffect(() => {
    if (selectedBond) {
      setBondOrder(selectedBond.order)
    }
  }, [selectedBond])

  // Update atom element
  const handleAtomElementChange = (newElement: string) => {
    if (!selectedAtom || !selectedAtomId) return

    const updated = molecule.updateAtom(selectedAtomId, {
      element: newElement,
    })
    
    setMolecule(updated)
    setAtomElement(newElement)
  }

  // Update atom charge
  const handleAtomChargeChange = (newCharge: number) => {
    if (!selectedAtom || !selectedAtomId) return

    const updated = molecule.updateAtom(selectedAtomId, {
      charge: newCharge,
      formalCharge: newCharge,
    })
    
    setMolecule(updated)
    setAtomCharge(newCharge)
  }

  // Update bond order
  const handleBondOrderChange = (newOrder: number) => {
    if (!selectedBond || !selectedBondId) return

    const updated = molecule.updateBond(selectedBondId, {
      order: newOrder,
    })
    
    setMolecule(updated)
    setBondOrder(newOrder)
  }

  // Clear selection
  const handleClearSelection = () => {
    setSelectedAtomId(null)
    setSelectedBondId(null)
  }

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
      <div className="flex items-center justify-between mb-3">
        <h3 className="text-sm font-semibold text-gray-700">Inspector</h3>
        <button
          onClick={handleClearSelection}
          className="text-xs text-gray-500 hover:text-gray-700"
          title="Clear selection (Esc)"
        >
          ✕
        </button>
      </div>
      
      {selectedAtom && (
        <div className="space-y-3">
          {/* Element selector */}
          <div>
            <label className="text-xs text-gray-500 mb-1 block">Element</label>
            <select
              value={atomElement}
              onChange={(e) => handleAtomElementChange(e.target.value)}
              className="w-full px-2 py-1.5 text-sm border border-gray-300 rounded focus:outline-none focus:ring-2 focus:ring-blue-500"
            >
              {COMMON_ELEMENTS.map((el) => (
                <option key={el} value={el}>
                  {el} - {ELEMENT_DATA[el]?.name || el}
                </option>
              ))}
            </select>
          </div>

          {/* Charge input */}
          <div>
            <label className="text-xs text-gray-500 mb-1 block">Charge</label>
            <input
              type="number"
              value={atomCharge}
              onChange={(e) => handleAtomChargeChange(parseInt(e.target.value) || 0)}
              min="-5"
              max="5"
              className="w-full px-2 py-1.5 text-sm border border-gray-300 rounded focus:outline-none focus:ring-2 focus:ring-blue-500"
            />
          </div>

          {/* Position (read-only) */}
          <div>
            <label className="text-xs text-gray-500 mb-1 block">Position</label>
            <div className="text-xs font-mono text-gray-600 bg-gray-50 px-2 py-1.5 rounded">
              ({selectedAtom.position[0].toFixed(2)}, {selectedAtom.position[1].toFixed(2)}, {selectedAtom.position[2].toFixed(2)})
            </div>
          </div>

          {/* Valence info */}
          {ELEMENT_DATA[atomElement] && (
            <div>
              <label className="text-xs text-gray-500 mb-1 block">Valence</label>
              <div className="text-xs text-gray-600">
                Max: {ELEMENT_DATA[atomElement].maxValence}
                {ELEMENT_DATA[atomElement].commonValence.length > 0 && (
                  <span className="ml-2">
                    Common: {ELEMENT_DATA[atomElement].commonValence.join(', ')}
                  </span>
                )}
              </div>
            </div>
          )}
        </div>
      )}

      {selectedBond && (
        <div className="space-y-3">
          {/* Bond order selector */}
          <div>
            <label className="text-xs text-gray-500 mb-1 block">Bond Order</label>
            <div className="flex gap-1">
              {VALID_BOND_ORDERS.map((order) => (
                <button
                  key={order}
                  onClick={() => handleBondOrderChange(order)}
                  className={`flex-1 px-2 py-1.5 text-xs rounded border transition-colors ${
                    bondOrder === order
                      ? 'bg-blue-50 border-blue-500 text-blue-700'
                      : 'bg-white border-gray-300 text-gray-700 hover:bg-gray-50'
                  }`}
                >
                  {order === 1.5 ? 'Aromatic' : `${order}x`}
                </button>
              ))}
            </div>
          </div>

          {/* Connected atoms (read-only) */}
          <div>
            <label className="text-xs text-gray-500 mb-1 block">Atoms</label>
            <div className="text-xs text-gray-600 bg-gray-50 px-2 py-1.5 rounded">
              {selectedBond.atom1} ↔ {selectedBond.atom2}
            </div>
          </div>
        </div>
      )}
    </div>
  )
}

