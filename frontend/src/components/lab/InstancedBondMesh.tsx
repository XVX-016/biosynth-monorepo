import React, { useRef, useMemo } from 'react'
import { useFrame } from '@react-three/fiber'
import * as THREE from 'three'
import type { Bond, Atom } from '../../types/molecule'

type Props = { bonds: Bond[], atoms: Atom[] }

export default function InstancedBondMesh({ bonds, atoms }: Props) {
  const ref = useRef<THREE.InstancedMesh | null>(null)
  const tmp = useMemo(() => new THREE.Object3D(), [])

  useFrame(() => {
    if (!ref.current) return
    bonds.forEach((b, i) => {
      const a1 = atoms.find(x => x.id === b.atom1)
      const a2 = atoms.find(x => x.id === b.atom2)
      if (!a1 || !a2) return
      
      const s = new THREE.Vector3(...a1.position)
      const e = new THREE.Vector3(...a2.position)
      const mid = s.clone().add(e).multiplyScalar(0.5)
      const dir = e.clone().sub(s)
      const len = dir.length()
      
      // Bond thickness based on order
      const radius = b.order === 1 ? 0.18 : b.order === 2 ? 0.22 : 0.26
      
      tmp.position.copy(mid)
      tmp.scale.set(radius, len, radius)
      tmp.lookAt(e)
      tmp.updateMatrix()
      ref.current!.setMatrixAt(i, tmp.matrix)
    })
    ref.current.count = bonds.length
    ref.current.instanceMatrix.needsUpdate = true
  })

  return (
    <instancedMesh ref={ref} args={[undefined, undefined, Math.max(2000, bonds.length + 100)]}>
      <cylinderGeometry args={[0.16, 0.16, 1, 12]} />
      <meshStandardMaterial color={0x9aa0a6} metalness={0.45} roughness={0.25} />
    </instancedMesh>
  )
}

