import React from 'react'
import * as THREE from 'three'

interface AtomMeshProps {
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

export default function AtomMesh({ position, element }: AtomMeshProps) {
  const radius = ELEMENT_RADII[element] || 1.0
  const color = ELEMENT_COLORS[element] || 0x9da3ae
  
  return (
    <mesh position={position}>
      <sphereGeometry args={[radius, 32, 32]} />
      <meshStandardMaterial 
        color={color} 
        metalness={0.9} 
        roughness={0.35} 
      />
    </mesh>
  )
}

