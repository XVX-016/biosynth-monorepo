// physicsWorker.ts - runs inside Worker context
type AtomState = { id: string; pos: [number, number, number] }
type BondState = { a: string; b: string }

self.onmessage = (ev: MessageEvent) => {
  const { cmd, atoms, bonds, iterations = 60 } = ev.data
  if (cmd === 'relax') {
    // simple spring-relax algorithm (very small)
    const positions: Record<string, number[]> = {}
    atoms.forEach((a: AtomState) => {
      positions[a.id] = [...a.pos]
    })
    
    const k = 0.12 // spring const
    const repulsion = 0.02
    
    for (let it = 0; it < iterations; it++) {
      // bond springs
      bonds.forEach((b: BondState) => {
        const p1 = positions[b.a]
        const p2 = positions[b.b]
        if (!p1 || !p2) return
        
        const dx = p2[0] - p1[0]
        const dy = p2[1] - p1[1]
        const dz = p2[2] - p1[2]
        const dist = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1e-6
        const desired = 1.5 // nominal bond length
        const diff = dist - desired
        const fx = (dx / dist) * diff * k
        const fy = (dy / dist) * diff * k
        const fz = (dz / dist) * diff * k
        
        positions[b.a][0] += fx
        positions[b.a][1] += fy
        positions[b.a][2] += fz
        positions[b.b][0] -= fx
        positions[b.b][1] -= fy
        positions[b.b][2] -= fz
      })
      
      // repulsion
      const ids = Object.keys(positions)
      for (let i = 0; i < ids.length; i++) {
        for (let j = i + 1; j < ids.length; j++) {
          const a = positions[ids[i]]
          const b = positions[ids[j]]
          if (!a || !b) continue
          
          const dx = a[0] - b[0]
          const dy = a[1] - b[1]
          const dz = a[2] - b[2]
          const dist = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1e-6
          const f = repulsion / (dist * dist)
          
          positions[ids[i]][0] += dx * f
          positions[ids[i]][1] += dy * f
          positions[ids[i]][2] += dz * f
          positions[ids[j]][0] -= dx * f
          positions[ids[j]][1] -= dy * f
          positions[ids[j]][2] -= dz * f
        }
      }
    }
    
    // return positions
    const out = Object.entries(positions).map(([id, pos]) => ({ id, pos }))
    self.postMessage({ cmd: 'relaxed', coords: out })
  }
}

