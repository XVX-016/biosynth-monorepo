import React from 'react'
import type { Atom } from '../../types/molecule'
import elementPalette from '../../data/elementPalette'
import { useLabStore } from '../../store/labStore'

interface AtomMeshProps {
  atoms: Atom[]
}

export default function AtomMesh({ atoms }: AtomMeshProps) {
  const selectedAtomId = useLabStore(s => s.selectedAtomId)
  
  // simple non-instanced version for foundation — later we'll swap to InstancedMesh
  return (
    <group>
      {atoms.map(a => {
        const palette = elementPalette[a.element] || elementPalette.C
        const isSelected = a.id === selectedAtomId
        return (
          <mesh 
            key={a.id} 
            position={a.position as any}
            userData={{ atomId: a.id }}
          >
            <sphereGeometry args={[palette.radius, 32, 24]} />
            <meshStandardMaterial
              color={isSelected ? palette.accent : palette.base}
              metalness={0.4}
              roughness={0.2}
            />
          </mesh>
        )
      })}
    </group>
  )
}

