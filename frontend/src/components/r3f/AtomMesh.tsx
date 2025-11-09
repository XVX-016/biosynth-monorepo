import React, { useRef, useState } from 'react'
import * as THREE from 'three'
import { useFrame } from '@react-three/fiber'
import { Outlines } from '@react-three/drei'
import { selectionManager } from './SelectionManager'
import { useMoleculeStore } from '../../store/moleculeStore'

interface AtomMeshProps {
  id: string
  position: [number, number, number]
  element: 'C' | 'H' | 'O' | 'N' | 'F' | 'S' | 'P' | 'Cl' | 'Br' | 'I'
}

const ELEMENT_RADII: Record<string, number> = {
  H: 0.55,
  C: 1.0,
  O: 0.9,
  N: 0.95,
  F: 0.9,
  S: 1.2,
  P: 1.15,
  Cl: 1.1,
  Br: 1.25,
  I: 1.4,
}

const ELEMENT_COLORS: Record<string, number> = {
  H: 0x9fb6ff,
  C: 0x9da3ae,
  O: 0xff6b6b,
  N: 0x8b5cf6,
  F: 0x10b981,
  S: 0xf59e0b,
  P: 0xf97316,
  Cl: 0x84cc16,
  Br: 0x991b1b,
  I: 0x7c3aed,
}

export default function AtomMesh({ id, position, element }: AtomMeshProps) {
  const meshRef = useRef<THREE.Mesh>(null)
  const [isHovered, setIsHovered] = useState(false)
  const [isDragging, setIsDragging] = useState(false)
  const selectedAtomId = useMoleculeStore((state) => state.selectedAtomId)
  const tool = useMoleculeStore((state) => state.tool)
  
  const radius = ELEMENT_RADII[element] || 1.0
  const color = ELEMENT_COLORS[element] || 0x9da3ae
  const isSelected = selectedAtomId === id
  
  // Visual feedback based on tool
  let scale = isSelected ? 1.15 : 1.0
  let outlineColor = 0x4EA7FF // Blue for select/hover
  if (tool === 'delete' && (isHovered || isSelected)) {
    outlineColor = 0xFF6B6B // Red for delete
  } else if (tool === 'bond' && (isHovered || isSelected)) {
    outlineColor = 0x4EA7FF // Blue for bond
  }
  
  const opacity = isDragging ? 0.7 : 1.0

  // Handle pointer over (hover)
  const handlePointerOver = (e: React.PointerEvent) => {
    e.stopPropagation()
    setIsHovered(true)
    selectionManager.onHover(id)
  }

  // Handle pointer out (unhover)
  const handlePointerOut = (e: React.PointerEvent) => {
    e.stopPropagation()
    setIsHovered(false)
    if (!isDragging) {
      selectionManager.onHover(null)
    }
  }

  // Handle click (select or delete)
  const handleClick = (e: React.MouseEvent) => {
    e.stopPropagation()
    
    if (tool === 'delete') {
      // Import removeAtom dynamically to avoid circular dependency
      import('../../lib/engineAdapter').then(({ removeAtom: remove }) => {
        remove(id)
      })
    } else {
      selectionManager.onSelect(id)
    }
  }

  // Handle pointer down (start drag)
  const handlePointerDown = (e: React.PointerEvent) => {
    e.stopPropagation()
    setIsDragging(true)
    selectionManager.startDrag(id)
    // Prevent orbit controls
    ;(e.target as HTMLElement).setPointerCapture(e.pointerId)
  }

  // Handle pointer move (position updates handled by InteractionLayer)
  const handlePointerMove = (e: React.PointerEvent) => {
    if (!isDragging) return
    e.stopPropagation()
    // Position updates are handled by InteractionLayer component
    // which has access to camera and can properly convert coordinates
  }

  // Handle pointer up (end drag)
  const handlePointerUp = (e: React.PointerEvent) => {
    if (isDragging) {
      e.stopPropagation()
      setIsDragging(false)
      selectionManager.endDrag()
      ;(e.target as HTMLElement).releasePointerCapture(e.pointerId)
    }
  }

  // Update position from store
  useFrame(() => {
    if (meshRef.current) {
      meshRef.current.position.set(...position)
    }
  })

  return (
    <mesh
      ref={meshRef}
      position={position}
      scale={scale}
      onClick={handleClick}
      onPointerOver={handlePointerOver}
      onPointerOut={handlePointerOut}
      onPointerDown={handlePointerDown}
      onPointerMove={handlePointerMove}
      onPointerUp={handlePointerUp}
      cursor={isDragging ? 'grabbing' : 'pointer'}
    >
      <sphereGeometry args={[radius, 64, 64]} />
      <meshPhysicalMaterial
        color={color}
        roughness={0.15}
        metalness={0.6}
        clearcoat={1}
        clearcoatRoughness={0.1}
        thickness={0.8}
        transparent={isDragging}
        opacity={opacity}
      />
      {(isHovered || isSelected) && (
        <Outlines thickness={0.1} color={outlineColor} />
      )}
    </mesh>
  )
}

