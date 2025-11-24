import React, { useEffect, useRef, useState } from 'react'
import { useLabStore } from '../../store/labStore'

export default function ExploreModeUI() {
  const mol = useLabStore(s => s.molecule)
  const [running, setRunning] = useState(false)
  const workerRef = useRef<Worker | null>(null)

  useEffect(() => {
    // Create worker - use inline blob for compatibility
    const workerCode = `
      self.onmessage = (ev) => {
        const { cmd, atoms, bonds, iterations = 60 } = ev.data
        if (cmd === 'relax') {
          const positions = {}
          atoms.forEach((a) => positions[a.id] = [...a.pos])
          const k = 0.12
          const repulsion = 0.02
          for (let it = 0; it < iterations; it++) {
            bonds.forEach((b) => {
              const p1 = positions[b.a], p2 = positions[b.b]
              if (!p1 || !p2) return
              const dx = p2[0] - p1[0], dy = p2[1] - p1[1], dz = p2[2] - p1[2]
              const dist = Math.sqrt(dx*dx+dy*dy+dz*dz) || 1e-6
              const desired = 1.5
              const diff = dist - desired
              const fx = (dx/dist) * diff * k, fy = (dy/dist) * diff * k, fz = (dz/dist) * diff * k
              positions[b.a][0] += fx; positions[b.a][1] += fy; positions[b.a][2] += fz
              positions[b.b][0] -= fx; positions[b.b][1] -= fy; positions[b.b][2] -= fz
            })
            const ids = Object.keys(positions)
            for (let i = 0; i < ids.length; i++) {
              for (let j = i + 1; j < ids.length; j++) {
                const a = positions[ids[i]], b = positions[ids[j]]
                if (!a || !b) continue
                const dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2]
                const dist = Math.sqrt(dx*dx+dy*dy+dz*dz) || 1e-6
                const f = repulsion/(dist*dist)
                positions[ids[i]][0] += dx*f; positions[ids[i]][1] += dy*f; positions[ids[i]][2] += dz*f
                positions[ids[j]][0] -= dx*f; positions[ids[j]][1] -= dy*f; positions[ids[j]][2] -= dz*f
              }
            }
          }
          const out = Object.entries(positions).map(([id,pos]) => ({ id, pos }))
          self.postMessage({ cmd:'relaxed', coords: out })
        }
      }
    `
    
    try {
      const blob = new Blob([workerCode], { type: 'application/javascript' })
      const worker = new Worker(URL.createObjectURL(blob))
      workerRef.current = worker
      
      worker.onmessage = (ev: MessageEvent) => {
        if (ev.data.cmd === 'relaxed') {
          // apply positions into store with animation
          ev.data.coords.forEach((c: any) => {
            useLabStore.getState().moveAtom(c.id, [c.pos[0], c.pos[1], c.pos[2]])
          })
          setRunning(false)
        }
      }
    } catch (e) {
      console.error('Failed to create worker:', e)
    }
    
    return () => {
      if (workerRef.current) {
        workerRef.current.terminate()
        workerRef.current = null
      }
    }
  }, [])

  const relaxClientSide = () => {
    if (!mol.atoms.length || !workerRef.current) return
    setRunning(true)
    const atoms = mol.atoms.map(a => ({ id: a.id, pos: a.position }))
    const bonds = mol.bonds.map(b => ({ a: b.atom1, b: b.atom2 }))
    workerRef.current.postMessage({ cmd: 'relax', atoms, bonds, iterations: 80 })
  }

  const relaxServer = async () => {
    // call backend /api/relax which will return optimized coords (later RDKit)
    try {
      const resp = await fetch('/api/relax', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(mol)
      })
      if (!resp.ok) {
        console.error('relax failed')
        return
      }
      const json = await resp.json()
      if (json.atoms) {
        json.atoms.forEach((a: any) => {
          useLabStore.getState().moveAtom(a.id, [a.x, a.y, a.z])
        })
      }
    } catch (error) {
      console.error('Server relax error:', error)
    }
  }

  return (
    <div className="flex flex-col gap-2">
      <button
        onClick={relaxClientSide}
        disabled={running || !workerRef.current}
        className="w-full px-2 py-1.5 text-xs rounded border border-gray-300 bg-white hover:bg-gray-50 disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
      >
        {running ? 'Relaxing...' : 'Quick Relax'}
      </button>
      <button 
        onClick={relaxServer} 
        className="w-full px-2 py-1.5 text-xs rounded border border-gray-300 bg-white hover:bg-gray-50 transition-colors"
      >
        Server Relax
      </button>
    </div>
  )
}
