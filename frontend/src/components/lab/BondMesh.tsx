import React from 'react'
import * as THREE from 'three'
import type { Bond, Atom } from '../../types/molecule'

interface BondMeshProps {
  bonds: Bond[]
  atoms: Atom[]
}

export default function BondMesh({ bonds, atoms }: BondMeshProps) {
  return (
    <group>
      {bonds.map(b => {
        const a1 = atoms.find(a => a.id === b.atom1)
        const a2 = atoms.find(a => a.id === b.atom2)
        if (!a1 || !a2) return null

        const start = new THREE.Vector3(...a1.position)
        const end = new THREE.Vector3(...a2.position)
        const dir = new THREE.Vector3().subVectors(end, start)
        const len = dir.length()
        const mid = start.clone().add(end).multiplyScalar(0.5)
        const quaternion = new THREE.Quaternion().setFromUnitVectors(
          new THREE.Vector3(0, 1, 0),
          dir.clone().normalize()
        )

        // Bond thickness based on order
        const radius = b.order === 1 ? 0.18 : b.order === 2 ? 0.22 : 0.26

        return (
          <mesh 
            key={b.id} 
            position={mid.toArray()} 
            quaternion={quaternion}
            userData={{ bondId: b.id }}
          >
            <cylinderGeometry args={[radius, radius, len, 12]} />
            <meshStandardMaterial color={0x9aa0a6} metalness={0.5} roughness={0.25} />
          </mesh>
        )
      })}
    </group>
  )
}

