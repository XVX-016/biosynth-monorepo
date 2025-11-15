import React, { useEffect, useMemo, useState } from 'react'
import { motion } from 'framer-motion'
import MoleculeViewer from '../components/MoleculeViewer'
import HeroScene from '../components/home/HeroScene'
import { useMoleculeStore } from '../store/moleculeStore'
import { listMolecules, MoleculeItem } from '../lib/api'
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
    <motion.div initial={{ opacity: 0 }} animate={{ opacity: 1 }} transition={{ duration: 0.3 }} className="space-y-6">
      {/* Hero Scene */}
      <div className="w-full h-[400px] rounded-xl overflow-hidden border border-chrome/20 shadow-glass">
        <HeroScene molecule={currentMolecule} />
      </div>
      
      <div className="grid grid-cols-12 gap-4">
        <div className="col-span-12 lg:col-span-8 space-y-4">
          <div className="grid grid-cols-1 sm:grid-cols-3 gap-4">
            <Card header={<div className="text-sm text-chrome">Total saved</div>}>
              <div className="text-3xl font-bold text-ivory">{kpis.total}</div>
            </Card>
            <Card header={<div className="text-sm text-chrome">Avg stability</div>}>
              <div className="text-3xl font-bold text-ivory">{kpis.avgStability}</div>
              <div className="mt-2 h-6">
                <svg viewBox="0 0 100 24" className="w-full h-full text-neonCyan">
                  <polyline fill="none" stroke="currentColor" strokeWidth="2" points="0,20 10,18 20,12 30,14 40,9 50,10 60,7 70,9 80,6 90,8 100,4" />
                </svg>
              </div>
            </Card>
            <Card header={<div className="text-sm text-chrome">Last generated</div>}>
              <div className="text-3xl font-bold truncate text-ivory" title={kpis.lastGenerated}>{kpis.lastGenerated}</div>
            </Card>
          </div>

          <Card header={<div className="flex items-center justify-between"><span className="font-semibold text-ivory">Recent Molecules</span></div>}>
            <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-4">
              {recent.slice(0, 6).map((m) => (
                <div key={m.id} className="rounded-xl border border-chrome/20 overflow-hidden frosted-glass hover:border-neonCyan/30 transition-colors">
                  <div className="h-28 flex items-center justify-center overflow-hidden bg-spaceGrey">
                    {m.thumbnail_b64 ? (
                      <img src={m.thumbnail_b64} alt={m.name} className="max-h-full max-w-full object-contain" />
                    ) : (
                      <div className="text-chrome text-sm">No preview</div>
                    )}
                  </div>
                  <div className="px-3 py-2">
                    <div className="font-medium truncate text-ivory">{m.name}</div>
                    <div className="text-xs text-chrome">{new Date(m.created_at).toLocaleDateString()}</div>
                  </div>
                </div>
              ))}
              {recent.length === 0 && <div className="text-chrome">No molecules yet.</div>}
            </div>
          </Card>
        </div>

        <div className="col-span-12 lg:col-span-4">
          <Card header={<div className="font-semibold text-ivory">3D Preview</div>}>
            <div className="h-[320px] rounded-lg border border-chrome/20 p-1 bg-spaceGrey">
              <MoleculeViewer />
            </div>
          </Card>
        </div>
      </div>
    </motion.div>
  )
}

