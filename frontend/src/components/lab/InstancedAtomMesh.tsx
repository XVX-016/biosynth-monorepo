import React, { useRef, useMemo } from 'react'
import { useFrame } from '@react-three/fiber'
import * as THREE from 'three'
import type { Atom } from '../../types/molecule'
import elementPalette from '../../data/elementPalette'
import { useLabStore } from '../../store/labStore'

type Props = { atoms: Atom[] }

export default function InstancedAtomMesh({ atoms }: Props) {
  const meshRef = useRef<THREE.InstancedMesh | null>(null)
  const tmpObj = useMemo(() => new THREE.Object3D(), [])
  const selectedAtomId = useLabStore(s => s.selectedAtomId)
  
  // CPK element colors for clean lab aesthetic
  const cpkColors: Record<string, number> = {
    H: 0xffffff,   // White
    C: 0x909090,   // Dark gray
    N: 0x3050f8,   // Blue
    O: 0xff0d0d,   // Red
    F: 0x90e050,   // Green
    P: 0xff8000,   // Orange
    S: 0xffff30,   // Yellow
    Cl: 0x1ff01f,  // Green
    Br: 0xa62929,  // Dark red
    I: 0x940094,   // Purple
  }

  useFrame(() => {
    if (!meshRef.current) return
    atoms.forEach((a, i) => {
      const palette = elementPalette[a.element] || elementPalette.C
      const radius = palette.radius || 1.0
      
      tmpObj.position.set(a.position[0], a.position[1], a.position[2])
      tmpObj.scale.setScalar(radius)
      tmpObj.updateMatrix()
      meshRef.current!.setMatrixAt(i, tmpObj.matrix)
      
      // Use CPK colors - flat, no metallic
      const isSelected = a.id === selectedAtomId
      const cpkColor = cpkColors[a.element] || 0xcccccc
      const color = isSelected 
        ? new THREE.Color(0x4676ff) // Blue highlight when selected
        : new THREE.Color(cpkColor)
      meshRef.current!.setColorAt(i, color)
    })
    meshRef.current.count = atoms.length
    meshRef.current.instanceMatrix.needsUpdate = true
    // @ts-ignore: instanceColor property exists
    if (meshRef.current.instanceColor) {
      meshRef.current.instanceColor.needsUpdate = true
    }
  })

  return (
    <instancedMesh ref={meshRef} args={[undefined, undefined, Math.max(2000, atoms.length + 100)]}>
      <sphereGeometry args={[1.0, 32, 20]} />
      <meshStandardMaterial metalness={0.0} roughness={0.6} />
    </instancedMesh>
  )
}
