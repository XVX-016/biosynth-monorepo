import React, { useEffect, useRef } from 'react'
import * as THREE from 'three'
import gsap from 'gsap'
import type { Bond, Atom } from '../../types/molecule'

interface AnimatedBondProps {
  bond: Bond
  atom1: Atom
  atom2: Atom
  index: number
}

export default function AnimatedBond({ bond, atom1, atom2, index }: AnimatedBondProps) {
  const meshRef = useRef<THREE.Mesh>(null)

  // Validate data
  if (!atom1 || !atom2 || !atom1.position || !atom2.position || 
      !Array.isArray(atom1.position) || !Array.isArray(atom2.position) ||
      atom1.position.length !== 3 || atom2.position.length !== 3) {
    return null
  }

  useEffect(() => {
    if (meshRef.current) {
      // Grow animation
      gsap.fromTo(
        meshRef.current.scale,
        { y: 0 },
        {
          y: 1,
          duration: 0.25,
          delay: index * 0.02,
          ease: 'power2.out',
        }
      )
    }
  }, [index])

  const start = new THREE.Vector3(
    Number(atom1.position[0]) || 0,
    Number(atom1.position[1]) || 0,
    Number(atom1.position[2]) || 0
  )
  const end = new THREE.Vector3(
    Number(atom2.position[0]) || 0,
    Number(atom2.position[1]) || 0,
    Number(atom2.position[2]) || 0
  )
  const dir = new THREE.Vector3().subVectors(end, start)
  const len = dir.length()
  
  if (len === 0) return null // Skip zero-length bonds
  
  const mid = start.clone().add(end).multiplyScalar(0.5)
  const quaternion = new THREE.Quaternion().setFromUnitVectors(
    new THREE.Vector3(0, 1, 0),
    dir.clone().normalize()
  )

  const radius = bond.order === 1 ? 0.18 : bond.order === 2 ? 0.22 : 0.26

  return (
    <mesh
      ref={meshRef}
      position={mid.toArray()}
      quaternion={quaternion}
      userData={{ bondId: bond.id }}
    >
      <cylinderGeometry args={[radius, radius, len, 12]} />
      <meshStandardMaterial color={0x9aa0a6} metalness={0.0} roughness={0.6} />
    </mesh>
  )
}

