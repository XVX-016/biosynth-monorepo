import React, { useEffect, useMemo, useState } from 'react'
import { motion } from 'framer-motion'
import MoleculeViewer from '../components/MoleculeViewer'
import HeroScene from '../components/home/HeroScene'
import { useMoleculeStore } from '../store/moleculeStore'
import { listMolecules } from '../lib/api'
import type { MoleculeItem } from '../lib/api'
import Card from '../components/ui/Card'

export default function Dashboard(){
  const currentMolecule = useMoleculeStore((state) => state.currentMolecule)
  const [recent, setRecent] = useState<MoleculeItem[]>([])

  useEffect(() => {
    let cancelled = false
    ;(async () => {
      try {
        const items = await listMolecules(12)
        if (!cancelled) setRecent(items ?? [])
      } catch {
        if (!cancelled) setRecent([])
      }
    })()
    return () => { cancelled = true }
  }, [])

  const kpis = useMemo(() => {
    const total = recent.length
    const avgStability = 0.86 // placeholder
    const lastGenerated = recent[0]?.name ?? '—'
    return { total, avgStability, lastGenerated }
  }, [recent])

  return (
    <motion.div initial={{ opacity: 0 }} animate={{ opacity: 1 }} transition={{ duration: 0.3 }} className="space-y-8">
      {/* Hero Section */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-8 items-center">
        <div className="space-y-6">
          <h1 className="text-5xl font-bold text-black">Design the Future of Molecular Science</h1>
          <p className="text-xl text-darkGrey leading-relaxed">
            Build, analyze, and explore molecular structures with cutting-edge AI-powered tools.
          </p>
          <div className="flex gap-4">
            <button className="btn-primary px-6 py-3">
              Generate Molecule
            </button>
            <button className="btn-secondary px-6 py-3">
              Explore Library
            </button>
          </div>
        </div>
        <div className="w-full h-[400px] rounded-xl overflow-hidden border border-lightGrey shadow-neon bg-white">
          <HeroScene molecule={currentMolecule} />
        </div>
      </div>
      
      {/* Stats Grid */}
      <div className="grid grid-cols-1 sm:grid-cols-3 gap-6">
        <Card className="p-6">
          <div className="text-sm text-midGrey mb-2">Total saved</div>
          <div className="text-4xl font-bold text-black">{kpis.total}</div>
        </Card>
        <Card className="p-6">
          <div className="text-sm text-midGrey mb-2">Avg stability</div>
          <div className="text-4xl font-bold text-black">{kpis.avgStability}</div>
          <div className="mt-3 h-6">
            <svg viewBox="0 0 100 24" className="w-full h-full text-darkGrey">
              <polyline fill="none" stroke="currentColor" strokeWidth="2" points="0,20 10,18 20,12 30,14 40,9 50,10 60,7 70,9 80,6 90,8 100,4" />
            </svg>
          </div>
        </Card>
        <Card className="p-6">
          <div className="text-sm text-midGrey mb-2">Last generated</div>
          <div className="text-2xl font-bold truncate text-black" title={kpis.lastGenerated}>{kpis.lastGenerated}</div>
        </Card>
      </div>

      {/* Recent Molecules */}
      <div>
        <h2 className="text-2xl font-semibold text-black mb-6">Recent Molecules</h2>
        <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-6">
          {recent.slice(0, 8).map((m) => (
            <Card key={m.id} className="overflow-hidden hover:shadow-neon-hover transition-shadow">
              <div className="h-32 flex items-center justify-center overflow-hidden bg-offwhite">
                {m.thumbnail_b64 ? (
                  <img src={m.thumbnail_b64} alt={m.name} className="max-h-full max-w-full object-contain" />
                ) : (
                  <div className="text-midGrey text-sm">No preview</div>
                )}
              </div>
              <div className="px-4 py-3">
                <div className="font-medium truncate text-black">{m.name}</div>
                <div className="text-xs text-midGrey mt-1">{new Date(m.created_at).toLocaleDateString()}</div>
              </div>
            </Card>
          ))}
          {recent.length === 0 && <div className="text-midGrey col-span-full text-center py-8">No molecules yet.</div>}
        </div>
      </div>

      {/* AI Model Cards */}
      <div>
        <h2 className="text-2xl font-semibold text-black mb-6">AI Models</h2>
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-6">
          <Card className="p-6 hover:shadow-neon-hover transition-shadow">
            <h3 className="text-lg font-semibold text-black mb-2">CoreGen-1</h3>
            <p className="text-sm text-darkGrey">Molecule generation from SMILES</p>
          </Card>
          <Card className="p-6 hover:shadow-neon-hover transition-shadow">
            <h3 className="text-lg font-semibold text-black mb-2">PropNet-X</h3>
            <p className="text-sm text-darkGrey">Property prediction</p>
          </Card>
          <Card className="p-6 hover:shadow-neon-hover transition-shadow">
            <h3 className="text-lg font-semibold text-black mb-2">ReactFlow-R1</h3>
            <p className="text-sm text-darkGrey">Synthesis planning</p>
          </Card>
          <Card className="p-6 hover:shadow-neon-hover transition-shadow">
            <h3 className="text-lg font-semibold text-black mb-2">ChemGPT-S</h3>
            <p className="text-sm text-darkGrey">Reaction simulation</p>
          </Card>
        </div>
      </div>
    </motion.div>
  )
}

