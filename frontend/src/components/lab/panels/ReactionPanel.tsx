import React, { useState } from 'react'
import { useLabStore } from '../../../store/labStore'
import { Line } from 'react-chartjs-2'
import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  PointElement,
  LineElement,
  Title,
  Tooltip,
  Legend
} from 'chart.js'

ChartJS.register(
  CategoryScale,
  LinearScale,
  PointElement,
  LineElement,
  Title,
  Tooltip,
  Legend
)

export default function ReactionPanel() {
  const molecule = useLabStore(s => s.molecule)
  const [reactionData, setReactionData] = useState<any>(null)
  const [mechanismData, setMechanismData] = useState<any>(null)
  const [loading, setLoading] = useState(false)
  const [activeTab, setActiveTab] = useState<'simulate' | 'mechanism'>('simulate')
  const [reactionType, setReactionType] = useState('substitution')

  const simulateReaction = async () => {
    if (molecule.atoms.length === 0) return
    
    setLoading(true)
    try {
      const res = await fetch('/api/reaction/simulate', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          reactants: [{
            atoms: molecule.atoms,
            bonds: molecule.bonds
          }],
          reaction_type: reactionType
        })
      })
      if (res.ok) {
        setReactionData(await res.json())
      }
    } catch (error) {
      console.error('Failed to simulate reaction:', error)
    } finally {
      setLoading(false)
    }
  }

  const predictMechanism = async () => {
    if (molecule.atoms.length === 0) return
    
    setLoading(true)
    try {
      const res = await fetch('/api/reaction/mechanism', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          reactants: [{
            atoms: molecule.atoms,
            bonds: molecule.bonds
          }],
          reaction_type: reactionType,
          max_steps: 5
        })
      })
      if (res.ok) {
        setMechanismData(await res.json())
      }
    } catch (error) {
      console.error('Failed to predict mechanism:', error)
    } finally {
      setLoading(false)
    }
  }

  const loadProduct = (product: any) => {
    if (product && product.atoms && product.bonds) {
      useLabStore.getState().loadMolecule(product)
    }
  }

  return (
    <div className="bg-white rounded-lg border border-gray-200 p-4 shadow-sm">
      <h4 className="text-sm font-semibold text-gray-700 mb-3">Reaction Simulation</h4>
      
      {molecule.atoms.length === 0 ? (
        <p className="text-xs text-gray-500">Add atoms to simulate reactions</p>
      ) : (
        <>
          {/* Reaction Type Selector */}
          <div className="mb-3">
            <label className="text-xs text-gray-600 mb-1 block">Reaction Type</label>
            <select
              value={reactionType}
              onChange={(e) => setReactionType(e.target.value)}
              className="w-full px-2 py-1 text-xs border border-gray-300 rounded"
            >
              <option value="substitution">Substitution</option>
              <option value="addition">Addition</option>
              <option value="elimination">Elimination</option>
            </select>
          </div>

          {/* Tabs */}
          <div className="flex gap-2 mb-4 border-b border-gray-200">
            <button
              onClick={() => setActiveTab('simulate')}
              className={`px-3 py-1 text-xs font-medium ${
                activeTab === 'simulate'
                  ? 'border-b-2 border-blue-600 text-blue-600'
                  : 'text-gray-600 hover:text-gray-900'
              }`}
            >
              Simulate
            </button>
            <button
              onClick={() => setActiveTab('mechanism')}
              className={`px-3 py-1 text-xs font-medium ${
                activeTab === 'mechanism'
                  ? 'border-b-2 border-blue-600 text-blue-600'
                  : 'text-gray-600 hover:text-gray-900'
              }`}
            >
              Mechanism
            </button>
          </div>

          {/* Simulate Tab */}
          {activeTab === 'simulate' && (
            <div className="space-y-3">
              <button
                onClick={simulateReaction}
                disabled={loading}
                className="w-full px-3 py-2 text-xs bg-blue-600 text-white rounded hover:bg-blue-700 disabled:opacity-50"
              >
                {loading ? 'Simulating...' : 'Simulate Reaction'}
              </button>
              
              {reactionData && (
                <div className="text-xs space-y-2">
                  <div className="flex justify-between">
                    <span className="text-gray-600">Feasible:</span>
                    <span className={reactionData.feasible ? 'text-green-600' : 'text-red-600'}>
                      {reactionData.feasible ? 'Yes' : 'No'}
                    </span>
                  </div>
                  
                  {reactionData.thermodynamics && (
                    <>
                      <div className="flex justify-between">
                        <span className="text-gray-600">ΔE:</span>
                        <span>{reactionData.thermodynamics.delta_E?.toFixed(2)} kcal/mol</span>
                      </div>
                      <div className="flex justify-between">
                        <span className="text-gray-600">ΔG:</span>
                        <span>{reactionData.thermodynamics.delta_G?.toFixed(2)} kcal/mol</span>
                      </div>
                    </>
                  )}
                  
                  {reactionData.products && reactionData.products.length > 0 && (
                    <div className="mt-3">
                      <div className="text-xs font-semibold text-gray-700 mb-2">Products:</div>
                      {reactionData.products.map((product: any, i: number) => (
                        <button
                          key={i}
                          onClick={() => loadProduct(product)}
                          className="w-full px-2 py-1 mb-1 text-xs bg-gray-100 hover:bg-gray-200 rounded text-left"
                        >
                          Product {i + 1} ({product.atoms?.length || 0} atoms)
                        </button>
                      ))}
                    </div>
                  )}
                </div>
              )}
            </div>
          )}

          {/* Mechanism Tab */}
          {activeTab === 'mechanism' && (
            <div className="space-y-3">
              <button
                onClick={predictMechanism}
                disabled={loading}
                className="w-full px-3 py-2 text-xs bg-blue-600 text-white rounded hover:bg-blue-700 disabled:opacity-50"
              >
                {loading ? 'Predicting...' : 'Predict Mechanism'}
              </button>
              
              {mechanismData && (
                <div className="text-xs space-y-2">
                  <div className="flex justify-between">
                    <span className="text-gray-600">Steps:</span>
                    <span>{mechanismData.steps?.length || 0}</span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-gray-600">Activation Energy:</span>
                    <span>{mechanismData.activation_energy?.toFixed(2)} kcal/mol</span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-gray-600">Reaction Energy:</span>
                    <span>{mechanismData.reaction_energy?.toFixed(2)} kcal/mol</span>
                  </div>
                  
                  {mechanismData.steps && mechanismData.steps.length > 0 && (
                    <>
                      <div className="mt-3 h-32">
                        <Line
                          data={{
                            labels: mechanismData.steps.map((s: any) => `Step ${s.step}`),
                            datasets: [{
                              label: 'Energy',
                              data: mechanismData.steps.map((s: any) => s.energy),
                              borderColor: '#4676ff',
                              tension: 0.4
                            }]
                          }}
                          options={{
                            responsive: true,
                            maintainAspectRatio: false,
                            scales: {
                              y: {
                                title: { display: true, text: 'Energy (kcal/mol)' }
                              },
                              x: {
                                title: { display: true, text: 'Step' }
                              }
                            }
                          }}
                        />
                      </div>
                      
                      <div className="mt-2 max-h-32 overflow-y-auto">
                        {mechanismData.steps.map((step: any, i: number) => (
                          <button
                            key={i}
                            onClick={() => loadProduct(step.intermediate)}
                            className="w-full px-2 py-1 mb-1 text-xs bg-gray-100 hover:bg-gray-200 rounded text-left"
                          >
                            {step.description} - {step.energy.toFixed(2)} kcal/mol
                          </button>
                        ))}
                      </div>
                    </>
                  )}
                </div>
              )}
            </div>
          )}
        </>
      )}
    </div>
  )
}

