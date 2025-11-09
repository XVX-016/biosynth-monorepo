import React, { useEffect, useState } from 'react'
import { motion } from 'framer-motion'
import MoleculeViewer from '../components/MoleculeViewer'
import ToolPanel from '../components/ToolPanel'
import AtomPalette from '../components/AtomPalette'
import { useMoleculeStore } from '../store/moleculeStore'
import { saveMolecule } from '../lib/api'
import { moleculeToJSON, getCanvasThumbnail } from '../lib/engineAdapter'
import { MoleculeSerializer } from '@biosynth/engine'

export default function Dashboard(){
  const generateMolecule = useMoleculeStore((state) => state.generateMolecule)
  const fetchPredictions = useMoleculeStore((state) => state.fetchPredictions)
  const currentMolecule = useMoleculeStore((state) => state.currentMolecule)
  const backendPredictions = useMoleculeStore((state) => state.backendPredictions)
  const loadingState = useMoleculeStore((state) => state.loadingState)
  const [saving, setSaving] = useState(false)

  // Auto-fetch predictions when molecule changes
  useEffect(() => {
    if (currentMolecule) {
      fetchPredictions()
    }
  }, [currentMolecule, fetchPredictions])

  const handleGenerate = async () => {
    await generateMolecule('Generate a new molecule')
  }

  const handleSave = async () => {
    if (!currentMolecule) {
      alert('No molecule to save')
      return
    }

    const name = prompt('Enter molecule name:')
    if (!name) return

    setSaving(true)
    try {
      const jsonGraph = moleculeToJSON(currentMolecule)
      const smiles = MoleculeSerializer.toSMILES(currentMolecule)
      const thumbnail = getCanvasThumbnail()
      const properties = backendPredictions ? JSON.stringify(backendPredictions) : null

      await saveMolecule({
        name,
        smiles: smiles || undefined,
        json_graph: jsonGraph,
        properties: properties || undefined,
        thumbnail_b64: thumbnail,
      })

      alert(`Molecule "${name}" saved successfully!`)
    } catch (error) {
      console.error('Failed to save molecule:', error)
      alert('Failed to save molecule')
    } finally {
      setSaving(false)
    }
  }

  return (
    <motion.div
      initial={{ opacity: 0 }}
      animate={{ opacity: 1 }}
      transition={{ duration: 0.3 }}
      className="space-y-6 relative"
    >
      <ToolPanel />
      <AtomPalette />
      <header className="flex items-center justify-between">
        <div className="text-xl font-semibold text-text-primary">BioSynth AI</div>
      </header>
      <section className="grid grid-cols-3 gap-6">
        <motion.div
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ delay: 0.1 }}
          className="col-span-2 p-8 bg-panel rounded-xl shadow-soft"
        >
          <h2 className="text-3xl font-bold">
            Design the Future of <span className="text-accent-blue">Molecular Science</span>
          </h2>
          <p className="text-text-secondary mt-2">
            Generate, analyze, and optimize molecular structures with cutting-edge AI.
          </p>
          <div className="mt-6">
            <div className="flex items-center justify-between mb-4">
              <h3 className="text-lg font-semibold text-text-primary">Molecule Lab</h3>
              <motion.button
                whileHover={{ scale: 1.05 }}
                whileTap={{ scale: 0.95 }}
                onClick={handleSave}
                disabled={!currentMolecule || saving}
                className="px-4 py-2 bg-accent-green text-white rounded-lg font-medium disabled:opacity-50 disabled:cursor-not-allowed"
              >
                {saving ? 'Saving...' : '💾 Save Molecule'}
              </motion.button>
            </div>
            <MoleculeViewer/>
          </div>
          {backendPredictions && (
            <div className="mt-4 p-4 bg-aluminum-light rounded-lg">
              <h3 className="font-semibold text-text-primary mb-2">Properties</h3>
              <div className="grid grid-cols-2 gap-2 text-sm">
                <div>Stability: {backendPredictions.stability?.toFixed(3)}</div>
                <div>Toxicity: {backendPredictions.toxicity?.toFixed(3)}</div>
                <div>Solubility: {backendPredictions.solubility?.toFixed(3)}</div>
                <div>Bioavailability: {backendPredictions.bioavailability?.toFixed(3)}</div>
                <div>Novelty: {backendPredictions.novelty?.toFixed(3)}</div>
              </div>
            </div>
          )}
          {loadingState === 'loading' && (
            <div className="mt-4 text-text-secondary">Loading predictions...</div>
          )}
        </motion.div>
        <motion.aside
          initial={{ opacity: 0, x: 20 }}
          animate={{ opacity: 1, x: 0 }}
          transition={{ delay: 0.2 }}
          className="p-8 bg-panel rounded-xl shadow-soft space-y-4"
        >
          <motion.button
            whileHover={{ scale: 1.05 }}
            whileTap={{ scale: 0.95 }}
            onClick={handleGenerate}
            disabled={loadingState === 'loading'}
            className="w-full bg-accent-blue text-white px-4 py-3 rounded-lg font-medium disabled:opacity-50"
          >
            {loadingState === 'loading' ? 'Generating...' : 'Generate Molecule'}
          </motion.button>
          <motion.button
            whileHover={{ scale: 1.05 }}
            whileTap={{ scale: 0.95 }}
            className="w-full bg-aluminum-DEFAULT px-4 py-3 rounded-lg font-medium text-text-primary"
          >
            Explore Library
          </motion.button>
        </motion.aside>
      </section>
    </motion.div>
  )
}

