import React, { useEffect, useState } from 'react'
import { motion } from 'framer-motion'
import { listMolecules, getMolecule, deleteMolecule, MoleculeItem } from '../lib/api'
import MoleculeCard from '../components/MoleculeCard'
import { useMoleculeStore } from '../store/moleculeStore'
import { moleculeFromJSON } from '../lib/engineAdapter'

export default function Library() {
  const [items, setItems] = useState<MoleculeItem[]>([])
  const [loading, setLoading] = useState(false)
  const setMolecule = useMoleculeStore((state) => state.setMolecule)

  useEffect(() => {
    loadMolecules()
  }, [])

  const loadMolecules = async () => {
    setLoading(true)
    try {
      const molecules = await listMolecules(50)
      setItems(molecules)
    } catch (error) {
      console.error('Failed to load molecules:', error)
    } finally {
      setLoading(false)
    }
  }

  const openInLab = async (id: number) => {
    try {
      const mol = await getMolecule(id)
      if (mol.json_graph) {
        const molecule = moleculeFromJSON(mol.json_graph)
        if (molecule) {
          setMolecule(molecule)
          // Navigate to dashboard (or lab page if you have one)
          window.location.href = '/'
        }
      }
    } catch (error) {
      console.error('Failed to load molecule:', error)
      alert('Failed to load molecule')
    }
  }

  const remove = async (id: number) => {
    if (!confirm('Delete this molecule?')) return
    
    try {
      await deleteMolecule(id)
      setItems(items.filter((i) => i.id !== id))
    } catch (error) {
      console.error('Failed to delete molecule:', error)
      alert('Failed to delete molecule')
    }
  }

  return (
    <motion.div
      initial={{ opacity: 0 }}
      animate={{ opacity: 1 }}
      transition={{ duration: 0.3 }}
      className="p-8 space-y-6"
    >
      <header className="flex items-center justify-between">
        <div>
          <h1 className="text-3xl font-bold text-text-primary">Molecule Library</h1>
          <p className="text-text-secondary mt-1">Your saved molecular structures</p>
        </div>
        <button
          onClick={loadMolecules}
          className="px-4 py-2 bg-accent-blue text-white rounded-lg font-medium hover:opacity-90"
        >
          Refresh
        </button>
      </header>

      {loading ? (
        <div className="text-center py-12 text-text-secondary">Loading molecules...</div>
      ) : items.length === 0 ? (
        <div className="text-center py-12">
          <div className="text-text-secondary mb-4">No molecules saved yet</div>
          <p className="text-text-tertiary">Create and save molecules in the Lab to see them here</p>
        </div>
      ) : (
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 xl:grid-cols-4 gap-6">
          {items.map((item) => (
            <MoleculeCard
              key={item.id}
              item={item}
              onOpen={() => openInLab(item.id)}
              onDelete={() => remove(item.id)}
            />
          ))}
        </div>
      )}
    </motion.div>
  )
}

