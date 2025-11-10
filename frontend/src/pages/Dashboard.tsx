import React, { useEffect, useMemo, useState } from 'react'
import { motion } from 'framer-motion'
import MoleculeViewer from '../components/MoleculeViewer'
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
      <div className="grid grid-cols-12 gap-4">
        <div className="col-span-12 lg:col-span-8 space-y-4">
          <div className="grid grid-cols-1 sm:grid-cols-3 gap-4">
            <Card header={<div className="text-sm text-text-secondary">Total saved</div>}>
              <div className="text-3xl font-bold">{kpis.total}</div>
            </Card>
            <Card header={<div className="text-sm text-text-secondary">Avg stability</div>}>
              <div className="text-3xl font-bold">{kpis.avgStability}</div>
              <div className="mt-2 h-6">
                <svg viewBox="0 0 100 24" className="w-full h-full text-accent-blue">
                  <polyline fill="none" stroke="currentColor" strokeWidth="2" points="0,20 10,18 20,12 30,14 40,9 50,10 60,7 70,9 80,6 90,8 100,4" />
                </svg>
              </div>
            </Card>
            <Card header={<div className="text-sm text-text-secondary">Last generated</div>}>
              <div className="text-3xl font-bold truncate" title={kpis.lastGenerated}>{kpis.lastGenerated}</div>
            </Card>
          </div>

          <Card header={<div className="flex items-center justify-between"><span className="font-semibold">Recent Molecules</span></div>}>
            <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-4">
              {recent.slice(0, 6).map((m) => (
                <div key={m.id} className="rounded-xl border border-aluminum-DEFAULT overflow-hidden bg-aluminum-light">
                  <div className="h-28 flex items-center justify-center overflow-hidden bg-white">
                    {m.thumbnail_b64 ? (
                      <img src={m.thumbnail_b64} alt={m.name} className="max-h-full max-w-full object-contain" />
                    ) : (
                      <div className="text-text-tertiary text-sm">No preview</div>
                    )}
                  </div>
                  <div className="px-3 py-2">
                    <div className="font-medium truncate">{m.name}</div>
                    <div className="text-xs text-text-tertiary">{new Date(m.created_at).toLocaleDateString()}</div>
                  </div>
                </div>
              ))}
              {recent.length === 0 && <div className="text-text-secondary">No molecules yet.</div>}
            </div>
          </Card>
        </div>

        <div className="col-span-12 lg:col-span-4">
          <Card header={<div className="font-semibold">3D Preview</div>}>
            <div className="h-[320px] rounded-lg border border-aluminum-DEFAULT p-1">
              <MoleculeViewer />
            </div>
          </Card>
        </div>
      </div>
    </motion.div>
  )
}

