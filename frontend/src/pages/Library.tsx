import { useEffect, useMemo, useState } from 'react'
import { motion } from 'framer-motion'
import { listMolecules, getMolecule, deleteMolecule } from '../lib/api'
import type { MoleculeItem } from '../lib/api'
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

  const { data: paged, totalPages } = useMemo(() => paginate(filtered, page, pageSize), [filtered, page, pageSize])
  
  // Only reset page when search query changes or when page exceeds totalPages after filtering
  useEffect(() => {
    const maxPage = Math.ceil(filtered.length / pageSize)
    if (maxPage > 0 && page > maxPage) {
      setPage(1)
    }
  }, [filtered.length, pageSize]) // Only depend on filtered length, not page or totalPages

  // Scroll to top when page changes
  useEffect(() => {
    window.scrollTo({ top: 0, behavior: 'smooth' })
  }, [page])

  return (
    <motion.div
      initial={{ opacity: 0 }}
      animate={{ opacity: 1 }}
      transition={{ duration: 0.3 }}
      className="p-8 space-y-6"
    >
      <motion.header 
        className="flex items-center justify-between gap-4 flex-wrap"
        initial={{ opacity: 0, y: -20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.4, delay: 0.1 }}
      >
        <div className="min-w-0">
          <h1 className="text-3xl font-bold text-black truncate">Molecule Library</h1>
          <p className="text-darkGrey mt-1">Your saved molecular structures</p>
        </div>
        <div className="flex items-center gap-3">
          <input
            value={q}
            onChange={(e) => { setQ(e.target.value); setPage(1); }}
            placeholder="Search by name or SMILES..."
            className="w-64 rounded-lg border border-lightGrey bg-white text-black px-3 py-2 outline-none focus:ring-2 focus:ring-darkGrey/20 focus:border-darkGrey placeholder:text-midGrey"
          />
          <button
            onClick={loadMolecules}
            className="btn-secondary px-4 py-2"
          >
            Refresh
          </button>
        </div>
      </motion.header>

      {loading ? (
        <div className="text-center py-12 text-midGrey">Loading molecules...</div>
      ) : filtered.length === 0 ? (
        <div className="text-center py-12">
          <div className="text-midGrey mb-2">No results</div>
          <p className="text-midGrey text-sm">Try clearing the search or saving molecules in the Lab</p>
        </div>
      ) : (
        <motion.div 
          key={`page-${page}`}
          className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6"
          initial={{ opacity: 0 }}
          animate={{ opacity: 1 }}
          transition={{ duration: 0.3 }}
        >
          {paged.map((item, index) => (
            <motion.div
              key={item.id}
              initial={{ opacity: 0, y: 20 }}
              animate={{ opacity: 1, y: 0 }}
              transition={{ duration: 0.2, delay: index * 0.03 }}
            >
              <MoleculeCard
                item={item}
                onOpen={() => openInLab(item.id)}
                onDelete={() => remove(item.id)}
              />
            </motion.div>
          ))}
        </motion.div>
      )}

      {/* Pagination */}
      {filtered.length > 0 && totalPages > 1 && (
        <motion.div 
          className="flex items-center justify-center gap-2 pt-2"
          initial={{ opacity: 0 }}
          animate={{ opacity: 1 }}
          transition={{ duration: 0.3, delay: 0.4 }}
        >
          <button
            onClick={() => setPage((p) => Math.max(1, p - 1))}
            disabled={page === 1}
            className="btn-secondary px-3 py-1 text-sm disabled:opacity-50"
          >
            Prev
          </button>
          <div className="text-sm text-darkGrey">
            Page <span className="font-semibold text-black">{page}</span> of {totalPages}
          </div>
          <button
            onClick={() => setPage((p) => Math.min(totalPages, p + 1))}
            disabled={page === totalPages}
            className="btn-secondary px-3 py-1 text-sm disabled:opacity-50"
          >
            Next
          </button>
        </motion.div>
      )}
    </motion.div>
  )
}

