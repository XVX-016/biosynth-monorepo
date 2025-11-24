import React, { useEffect, useRef, useState } from 'react'
import { useFrame } from '@react-three/fiber'
import * as THREE from 'three'
import gsap from 'gsap'
import type { Atom } from '../../types/molecule'
import { useLabStore } from '../../store/labStore'
import AtomTooltip from './AtomTooltip'

interface AnimatedAtomProps {
  atom: Atom
  index: number
}

// CPK colors
const cpkColors: Record<string, number> = {
  H: 0xffffff,
  C: 0x909090,
  N: 0x3050f8,
  O: 0xff0d0d,
  F: 0x90e050,
  P: 0xff8000,
  S: 0xffff30,
  Cl: 0x1ff01f,
  Br: 0xa62929,
  I: 0x940094,
}

export default function AnimatedAtom({ atom, index }: AnimatedAtomProps) {
  const meshRef = useRef<THREE.Mesh>(null)
  const [hovered, setHovered] = useState(false)
  const selectedAtomId = useLabStore(s => s.selectedAtomId)
  const isSelected = atom?.id === selectedAtomId

  // Validate atom data
  if (!atom || !atom.position || !Array.isArray(atom.position) || atom.position.length !== 3) {
    return null
  }

  // Ensure position is valid numbers
  const position: [number, number, number] = [
    Number(atom.position[0]) || 0,
    Number(atom.position[1]) || 0,
    Number(atom.position[2]) || 0,
  ]

  // Pop animation on mount
  useEffect(() => {
    if (meshRef.current) {
      gsap.fromTo(
        meshRef.current.scale,
        { x: 0, y: 0, z: 0 },
        {
          x: 1,
          y: 1,
          z: 1,
          duration: 0.4,
          delay: index * 0.05,
          ease: 'back.out(1.7)',
        }
      )
    }
  }, [index])

  // Selection highlight
  useEffect(() => {
    if (meshRef.current && isSelected) {
      gsap.to(meshRef.current.scale, {
        x: 1.2,
        y: 1.2,
        z: 1.2,
        duration: 0.2,
        ease: 'power2.out',
      })
    } else if (meshRef.current) {
      gsap.to(meshRef.current.scale, {
        x: 1,
        y: 1,
        z: 1,
        duration: 0.2,
        ease: 'power2.out',
      })
    }
  }, [isSelected])

  const color = isSelected ? 0x4676ff : (cpkColors[atom.element] || 0xcccccc)

  return (
    <group>
      <mesh
        ref={meshRef}
        position={position}
        userData={{ atomId: atom.id, type: 'atom' }}
        onPointerOver={(e) => {
          e.stopPropagation()
          setHovered(true)
        }}
        onPointerOut={() => setHovered(false)}
      >
        <sphereGeometry args={[1.0, 32, 20]} />
        <meshStandardMaterial
          color={color}
          metalness={0.0}
          roughness={0.6}
        />
      </mesh>
      {hovered && (
        <AtomTooltip atom={atom} position={position} />
      )}
    </group>
  )
}

