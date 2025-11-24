import React, { useState } from 'react'
import { useLabStore } from '../../../store/labStore'
import PathwayGraph from '../../retrosynthesis/PathwayGraph'
import PathwayExport from '../../retrosynthesis/PathwayExport'

export default function RetrosynthesisPanel() {
  const molecule = useLabStore(s => s.molecule)
  const [pathways, setPathways] = useState<any[]>([])
  const [loading, setLoading] = useState(false)
  const [selectedPathway, setSelectedPathway] = useState<number | null>(null)
  const [maxSteps, setMaxSteps] = useState(5)

  const planRetrosynthesis = async () => {
    if (molecule.atoms.length === 0) return
    
    setLoading(true)
    try {
      const res = await fetch('/api/retrosynthesis/plan', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          atoms: molecule.atoms,
          bonds: molecule.bonds,
          max_steps: maxSteps,
          max_pathways: 10
        })
      })
      if (res.ok) {
        const data = await res.json()
        setPathways(data.pathways || [])
        if (data.pathways && data.pathways.length > 0) {
          setSelectedPathway(0)
        }
      }
    } catch (error) {
      console.error('Failed to plan retrosynthesis:', error)
    } finally {
      setLoading(false)
    }
  }

  const loadStepMolecule = (step: any) => {
    if (step.molecule && step.molecule.atoms) {
      useLabStore.getState().loadMolecule(step.molecule)
    }
  }

  return (
    <div className="bg-white rounded-lg border border-gray-200 p-4 shadow-sm">
      <h4 className="text-sm font-semibold text-gray-700 mb-3">Retrosynthesis</h4>
      
      {molecule.atoms.length === 0 ? (
        <p className="text-xs text-gray-500">Add atoms to plan retrosynthesis</p>
      ) : (
        <>
          <div className="mb-3 space-y-2">
            <div>
              <label className="text-xs text-gray-600 mb-1 block">Max Steps</label>
              <input
                type="number"
                value={maxSteps}
                onChange={(e) => setMaxSteps(parseInt(e.target.value) || 5)}
                min={1}
                max={10}
                className="w-full px-2 py-1 text-xs border border-gray-300 rounded"
              />
            </div>
            
            <button
              onClick={planRetrosynthesis}
              disabled={loading}
              className="w-full px-3 py-2 text-xs bg-blue-600 text-white rounded hover:bg-blue-700 disabled:opacity-50"
            >
              {loading ? 'Planning...' : 'Plan Retrosynthesis'}
            </button>
          </div>

          {pathways.length > 0 && (
            <div className="space-y-3">
              <div className="text-xs font-semibold text-gray-700">
                Found {pathways.length} pathway{pathways.length !== 1 ? 's' : ''}
              </div>
              
              {/* Pathway Selector */}
              <div className="max-h-32 overflow-y-auto space-y-1">
                {pathways.map((pathway, index) => (
                  <button
                    key={index}
                    onClick={() => setSelectedPathway(index)}
                    className={`w-full px-2 py-1 text-xs rounded text-left ${
                      selectedPathway === index
                        ? 'bg-blue-100 border border-blue-300'
                        : 'bg-gray-100 hover:bg-gray-200'
                    }`}
                  >
                    Pathway {index + 1} - Score: {pathway.score?.toFixed(2) || 'N/A'} 
                    ({pathway.total_steps || 0} steps)
                  </button>
                ))}
              </div>

              {/* Selected Pathway Details */}
              {selectedPathway !== null && pathways[selectedPathway] && (
                <div className="mt-3 border-t border-gray-200 pt-3">
                  <div className="mb-2">
                    <PathwayExport pathway={pathways[selectedPathway]} />
                  </div>
                  <PathwayGraph
                    pathway={pathways[selectedPathway]}
                    onStepClick={loadStepMolecule}
                  />
                </div>
              )}
            </div>
          )}
        </>
      )}
    </div>
  )
}

