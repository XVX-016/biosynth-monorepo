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

export default function EnergyPanel() {
  const molecule = useLabStore(s => s.molecule)
  const [energyData, setEnergyData] = useState<any>(null)
  const [optimizationData, setOptimizationData] = useState<any>(null)
  const [simulationData, setSimulationData] = useState<any>(null)
  const [loading, setLoading] = useState(false)
  const [activeTab, setActiveTab] = useState<'calculate' | 'optimize' | 'simulate'>('calculate')

  const calculateEnergy = async () => {
    if (molecule.atoms.length === 0) return
    
    setLoading(true)
    try {
      const res = await fetch('/api/energy/calculate', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          atoms: molecule.atoms,
          bonds: molecule.bonds
        })
      })
      if (res.ok) {
        setEnergyData(await res.json())
      }
    } catch (error) {
      console.error('Failed to calculate energy:', error)
    } finally {
      setLoading(false)
    }
  }

  const optimizeGeometry = async () => {
    if (molecule.atoms.length === 0) return
    
    setLoading(true)
    try {
      const res = await fetch('/api/energy/optimize', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          atoms: molecule.atoms,
          bonds: molecule.bonds,
          max_iterations: 100
        })
      })
      if (res.ok) {
        const data = await res.json()
        setOptimizationData(data)
        
        // Update molecule with optimized coordinates
        if (data.optimized_molecule) {
          useLabStore.getState().loadMolecule(data.optimized_molecule)
        }
      }
    } catch (error) {
      console.error('Failed to optimize geometry:', error)
    } finally {
      setLoading(false)
    }
  }

  const runSimulation = async () => {
    if (molecule.atoms.length === 0) return
    
    setLoading(true)
    try {
      const res = await fetch('/api/energy/simulate', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          atoms: molecule.atoms,
          bonds: molecule.bonds,
          temperature: 300.0,
          steps: 100
        })
      })
      if (res.ok) {
        setSimulationData(await res.json())
      }
    } catch (error) {
      console.error('Failed to run simulation:', error)
    } finally {
      setLoading(false)
    }
  }

  return (
    <div className="bg-white rounded-lg border border-gray-200 p-4 shadow-sm">
      <h4 className="text-sm font-semibold text-gray-700 mb-3">Energy & Dynamics</h4>
      
      {molecule.atoms.length === 0 ? (
        <p className="text-xs text-gray-500">Add atoms to calculate energy</p>
      ) : (
        <>
          {/* Tabs */}
          <div className="flex gap-2 mb-4 border-b border-gray-200">
            <button
              onClick={() => setActiveTab('calculate')}
              className={`px-3 py-1 text-xs font-medium ${
                activeTab === 'calculate'
                  ? 'border-b-2 border-blue-600 text-blue-600'
                  : 'text-gray-600 hover:text-gray-900'
              }`}
            >
              Energy
            </button>
            <button
              onClick={() => setActiveTab('optimize')}
              className={`px-3 py-1 text-xs font-medium ${
                activeTab === 'optimize'
                  ? 'border-b-2 border-blue-600 text-blue-600'
                  : 'text-gray-600 hover:text-gray-900'
              }`}
            >
              Optimize
            </button>
            <button
              onClick={() => setActiveTab('simulate')}
              className={`px-3 py-1 text-xs font-medium ${
                activeTab === 'simulate'
                  ? 'border-b-2 border-blue-600 text-blue-600'
                  : 'text-gray-600 hover:text-gray-900'
              }`}
            >
              MD
            </button>
          </div>

          {/* Calculate Energy Tab */}
          {activeTab === 'calculate' && (
            <div className="space-y-3">
              <button
                onClick={calculateEnergy}
                disabled={loading}
                className="w-full px-3 py-2 text-xs bg-blue-600 text-white rounded hover:bg-blue-700 disabled:opacity-50"
              >
                {loading ? 'Calculating...' : 'Calculate Energy'}
              </button>
              
              {energyData && (
                <div className="text-xs space-y-1">
                  <div className="flex justify-between">
                    <span className="text-gray-600">Total Energy:</span>
                    <span className="font-semibold">{energyData.total_energy?.toFixed(2)} kcal/mol</span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-gray-600">Bond Energy:</span>
                    <span>{energyData.bond_energy?.toFixed(2)} kcal/mol</span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-gray-600">Angle Energy:</span>
                    <span>{energyData.angle_energy?.toFixed(2)} kcal/mol</span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-gray-600">Nonbonded:</span>
                    <span>{energyData.nonbonded_energy?.toFixed(2)} kcal/mol</span>
                  </div>
                </div>
              )}
            </div>
          )}

          {/* Optimize Geometry Tab */}
          {activeTab === 'optimize' && (
            <div className="space-y-3">
              <button
                onClick={optimizeGeometry}
                disabled={loading}
                className="w-full px-3 py-2 text-xs bg-blue-600 text-white rounded hover:bg-blue-700 disabled:opacity-50"
              >
                {loading ? 'Optimizing...' : 'Optimize Geometry'}
              </button>
              
              {optimizationData && (
                <div className="text-xs space-y-2">
                  <div className="flex justify-between">
                    <span className="text-gray-600">Initial Energy:</span>
                    <span>{optimizationData.initial_energy?.toFixed(2)} kcal/mol</span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-gray-600">Final Energy:</span>
                    <span className="font-semibold">{optimizationData.final_energy?.toFixed(2)} kcal/mol</span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-gray-600">Iterations:</span>
                    <span>{optimizationData.iterations}</span>
                  </div>
                  
                  {optimizationData.energy_history && optimizationData.energy_history.length > 1 && (
                    <div className="mt-3 h-32">
                      <Line
                        data={{
                          labels: optimizationData.energy_history.map((_: any, i: number) => i),
                          datasets: [{
                            label: 'Energy',
                            data: optimizationData.energy_history,
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
                              title: { display: true, text: 'Iteration' }
                            }
                          }
                        }}
                      />
                    </div>
                  )}
                </div>
              )}
            </div>
          )}

          {/* MD Simulation Tab */}
          {activeTab === 'simulate' && (
            <div className="space-y-3">
              <button
                onClick={runSimulation}
                disabled={loading}
                className="w-full px-3 py-2 text-xs bg-blue-600 text-white rounded hover:bg-blue-700 disabled:opacity-50"
              >
                {loading ? 'Simulating...' : 'Run MD Simulation'}
              </button>
              
              {simulationData && (
                <div className="text-xs space-y-2">
                  <div className="flex justify-between">
                    <span className="text-gray-600">Trajectory Steps:</span>
                    <span>{simulationData.trajectory?.length || 0}</span>
                  </div>
                  
                  {simulationData.energies && simulationData.energies.length > 0 && (
                    <div className="mt-3 h-32">
                      <Line
                        data={{
                          labels: simulationData.energies.map((e: any) => e.time.toFixed(2)),
                          datasets: [{
                            label: 'Energy',
                            data: simulationData.energies.map((e: any) => e.energy),
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
                              title: { display: true, text: 'Time (ps)' }
                            }
                          }
                        }}
                      />
                    </div>
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

