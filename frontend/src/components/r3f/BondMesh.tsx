import React, { useState } from 'react'
import * as THREE from 'three'
import { useMoleculeStore } from '../../store/moleculeStore'
import { Outlines } from '@react-three/drei'
import type { RenderMode } from '../../store/moleculeStore'

interface BondMeshProps {
  id: string
  from: [number, number, number]
  to: [number, number, number]
  order: number
  renderMode: RenderMode
  highlighted?: boolean
}

export default function BondMesh({ id, from, to, order, renderMode, highlighted = false }: BondMeshProps) {
  const [isHovered, setIsHovered] = useState(false)
  const selectedBondId = useMoleculeStore((state) => state.selectedBondId)
  const tool = useMoleculeStore((state) => state.tool)
  const selectBond = useMoleculeStore((state) => state.selectBond)
  
  const vFrom = new THREE.Vector3(...from)
  const vTo = new THREE.Vector3(...to)
  const diff = new THREE.Vector3().subVectors(vTo, vFrom)
  const length = diff.length()
  const mid = new THREE.Vector3().addVectors(vFrom, vTo).multiplyScalar(0.5)
  const q = new THREE.Quaternion().setFromUnitVectors(
    new THREE.Vector3(0, 1, 0), 
    diff.clone().normalize()
  )
  
  const isSelected = selectedBondId === id
  if (renderMode === 'spacefill') {
    return null
  }

  const baseRadius = order === 1 ? 0.14 : order === 2 ? 0.18 : 0.22
  const radius = renderMode === 'wireframe' ? baseRadius * 0.5 : baseRadius
  let outlineColor = 0x8BF3FF // neonCyan for select/hover
  
  if (tool === 'delete' && (isHovered || isSelected)) {
    outlineColor = 0xC6BDFE // violetEdge for delete
  }
  
  const handleClick = (e: React.MouseEvent) => {
    e.stopPropagation()
    
    if (tool === 'delete') {
      import('../../lib/engineAdapter').then(({ removeBond: remove }) => {
        remove(id)
      })
    } else {
      selectBond(id)
    }
  }
  
  return (
    <mesh 
      position={[mid.x, mid.y, mid.z]} 
      quaternion={[q.x, q.y, q.z, q.w]}
      onPointerOver={() => setIsHovered(true)}
      onPointerOut={() => setIsHovered(false)}
      onClick={handleClick}
      cursor={tool === 'delete' ? 'not-allowed' : 'pointer'}
    >
      <cylinderGeometry args={[radius, radius, length, 32]} />
      {renderMode === 'wireframe' ? (
        <meshBasicMaterial color={highlighted ? 0xffc857 : 0xb0b5c4} wireframe />
      ) : (
        <meshPhysicalMaterial
          color={
            highlighted
              ? 0xffc857
              : order === 1
              ? 0xc0c5d2
              : order === 2
              ? 0xd4d8e3
              : 0xe3c8c8
          }
          roughness={0.1}
          metalness={0.8}
          clearcoat={1}
          clearcoatRoughness={0.05}
          envMapIntensity={1.5}
        />
      )}
      {renderMode !== 'wireframe' && (isHovered || isSelected || highlighted) && (
        <Outlines thickness={0.05} color={outlineColor} />
      )}
    </mesh>
  )
}

