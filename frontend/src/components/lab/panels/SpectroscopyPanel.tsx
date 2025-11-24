import React, { useState, useEffect } from 'react'
import { useLabStore } from '../../../store/labStore'
import IRGraph from '../../spectroscopy/IRGraph'
import NMRGraph from '../../spectroscopy/NMRGraph'
import MassGraph from '../../spectroscopy/MassGraph'

export default function SpectroscopyPanel() {
  const molecule = useLabStore(s => s.molecule)
  const [irData, setIrData] = useState<any>(null)
  const [nmrData, setNmrData] = useState<any>(null)
  const [massData, setMassData] = useState<any>(null)
  const [loading, setLoading] = useState(false)
  const [activeTab, setActiveTab] = useState<'ir' | 'nmr' | 'mass'>('ir')
  const [highlightedAtoms, setHighlightedAtoms] = useState<string[]>([])

  useEffect(() => {
    if (molecule.atoms.length === 0) {
      setIrData(null)
      setNmrData(null)
      setMassData(null)
      return
    }

    const fetchSpectra = async () => {
      setLoading(true)
      try {
        // Fetch IR
        const irRes = await fetch('/api/spectroscopy/ir', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
            atoms: molecule.atoms,
            bonds: molecule.bonds
          })
        })
        if (irRes.ok) {
          setIrData(await irRes.json())
        }

        // Fetch NMR
        const nmrRes = await fetch('/api/spectroscopy/nmr', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
            atoms: molecule.atoms,
            bonds: molecule.bonds
          })
        })
        if (nmrRes.ok) {
          setNmrData(await nmrRes.json())
        }

        // Fetch Mass
        const massRes = await fetch('/api/spectroscopy/mass', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
            atoms: molecule.atoms,
            bonds: molecule.bonds
          })
        })
        if (massRes.ok) {
          setMassData(await massRes.json())
        }
      } catch (error) {
        console.error('Failed to fetch spectra:', error)
      } finally {
        setLoading(false)
      }
    }

    fetchSpectra()
  }, [molecule])

  const handlePeakHover = (peak: any) => {
    if (peak && peak.atoms) {
      setHighlightedAtoms(Array.isArray(peak.atoms) ? peak.atoms : [peak.atom])
    } else {
      setHighlightedAtoms([])
    }
  }

  return (
    <div className="bg-white rounded-lg border border-gray-200 p-4 shadow-sm">
      <h4 className="text-sm font-semibold text-gray-700 mb-3">Spectroscopy</h4>
      
      {molecule.atoms.length === 0 ? (
        <p className="text-xs text-gray-500">Add atoms to predict spectra</p>
      ) : loading ? (
        <p className="text-xs text-gray-500">Calculating spectra...</p>
      ) : (
        <>
          {/* Tabs */}
          <div className="flex gap-2 mb-4 border-b border-gray-200">
            <button
              onClick={() => setActiveTab('ir')}
              className={`px-3 py-1 text-xs font-medium ${
                activeTab === 'ir'
                  ? 'border-b-2 border-blue-600 text-blue-600'
                  : 'text-gray-600 hover:text-gray-900'
              }`}
            >
              IR
            </button>
            <button
              onClick={() => setActiveTab('nmr')}
              className={`px-3 py-1 text-xs font-medium ${
                activeTab === 'nmr'
                  ? 'border-b-2 border-blue-600 text-blue-600'
                  : 'text-gray-600 hover:text-gray-900'
              }`}
            >
              NMR
            </button>
            <button
              onClick={() => setActiveTab('mass')}
              className={`px-3 py-1 text-xs font-medium ${
                activeTab === 'mass'
                  ? 'border-b-2 border-blue-600 text-blue-600'
                  : 'text-gray-600 hover:text-gray-900'
              }`}
            >
              Mass
            </button>
          </div>

          {/* IR Tab */}
          {activeTab === 'ir' && irData && (
            <div>
              <IRGraph
                peaks={irData.peaks || []}
                onPeakHover={handlePeakHover}
              />
            </div>
          )}

          {/* NMR Tab */}
          {activeTab === 'nmr' && nmrData && (
            <div className="space-y-4">
              <div>
                <h5 className="text-xs font-semibold text-gray-700 mb-2">¹H NMR</h5>
                <NMRGraph
                  peaks={nmrData.protons || []}
                  type="1H"
                  onPeakHover={handlePeakHover}
                />
              </div>
              <div>
                <h5 className="text-xs font-semibold text-gray-700 mb-2">¹³C NMR</h5>
                <NMRGraph
                  peaks={nmrData.carbons || []}
                  type="13C"
                  onPeakHover={handlePeakHover}
                />
              </div>
            </div>
          )}

          {/* Mass Tab */}
          {activeTab === 'mass' && massData && (
            <div>
              <MassGraph
                peaks={massData.peaks || []}
                molecularIon={massData.molecular_ion || 0}
                onPeakHover={handlePeakHover}
              />
            </div>
          )}
        </>
      )}
    </div>
  )
}

