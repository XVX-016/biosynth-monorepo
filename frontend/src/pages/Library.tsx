import React, { useEffect, useMemo, useState } from 'react'
import { motion } from 'framer-motion'
import { listMolecules, getMolecule, deleteMolecule, MoleculeItem } from '../lib/api'
import MoleculeCard from '../components/MoleculeCard'
import { useMoleculeStore } from '../store/moleculeStore'
import { moleculeFromJSON } from '../lib/engineAdapter'
import { paginate } from '../lib/pagination'

export default function Library() {
  const [items, setItems] = useState<MoleculeItem[]>([])
  const [loading, setLoading] = useState(false)
  const [q, setQ] = useState('')
  const [page, setPage] = useState(1)
  const pageSize = 12
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

  const filtered = useMemo(() => {
    const query = q.trim().toLowerCase()
    if (!query) return items
    return items.filter(i => i.name.toLowerCase().includes(query) || (i.smiles || '').toLowerCase().includes(query))
  }, [items, q])

  const { data: paged, totalPages } = useMemo(() => paginate(filtered, page, pageSize), [filtered, page])
  useEffect(() => {
    // clamp page when filter changes
    if (page > totalPages) setPage(1)
  }, [totalPages])

  return (
    <motion.div
      initial={{ opacity: 0 }}
      animate={{ opacity: 1 }}
      transition={{ duration: 0.3 }}
      className="p-8 space-y-6"
    >
      <header className="flex items-center justify-between gap-4 flex-wrap">
        <div className="min-w-0">
          <h1 className="text-3xl font-bold text-ivory truncate">Molecule Library</h1>
          <p className="text-chrome mt-1">Your saved molecular structures</p>
        </div>
        <div className="flex items-center gap-3">
          <input
            value={q}
            onChange={(e) => { setQ(e.target.value); setPage(1); }}
            placeholder="Search by name or SMILES..."
            className="w-64 rounded-lg border border-chrome/20 bg-frostedGlass text-ivory px-3 py-2 outline-none focus:ring-2 focus:ring-neonCyan/50 focus:border-neonCyan/50 placeholder:text-chrome"
          />
          <button
            onClick={loadMolecules}
            className="px-4 py-2 bg-plasma-neon text-ionBlack rounded-lg font-medium hover:shadow-neon-sm transition-all"
          >
            Refresh
          </button>
        </div>
      </header>

      {loading ? (
        <div className="text-center py-12 text-chrome">Loading molecules...</div>
      ) : filtered.length === 0 ? (
        <div className="text-center py-12">
          <div className="text-chrome mb-2">No results</div>
          <p className="text-chrome/70 text-sm">Try clearing the search or saving molecules in the Lab</p>
        </div>
      ) : (
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 xl:grid-cols-4 gap-6">
          {paged.map((item) => (
            <MoleculeCard
              key={item.id}
              item={item}
              onOpen={() => openInLab(item.id)}
              onDelete={() => remove(item.id)}
            />
          ))}
        </div>
      )}

      {/* Pagination */}
      {filtered.length > 0 && totalPages > 1 && (
        <div className="flex items-center justify-center gap-2 pt-2">
          <button
            onClick={() => setPage((p) => Math.max(1, p - 1))}
            disabled={page === 1}
            className="px-3 py-1 rounded border border-chrome/20 bg-frostedGlass text-chrome hover:text-ivory hover:border-neonCyan/30 disabled:opacity-50 transition-all"
          >
            Prev
          </button>
          <div className="text-sm text-chrome">
            Page <span className="font-semibold text-ivory">{page}</span> of {totalPages}
          </div>
          <button
            onClick={() => setPage((p) => Math.min(totalPages, p + 1))}
            disabled={page === totalPages}
            className="px-3 py-1 rounded border border-chrome/20 bg-frostedGlass text-chrome hover:text-ivory hover:border-neonCyan/30 disabled:opacity-50 transition-all"
          >
            Next
          </button>
        </div>
      )}
    </motion.div>
  )
}

