/**
 * ThreeDViewer - Clean 3D molecule viewer
 * 
 * Phase 9: 3D Preview Integration
 * 
 * Features:
 * - Orbit controls
 * - Atom spheres (CPK colors)
 * - Bond cylinders
 * - Hover highlighting
 * - Lazy mount to avoid WebGL loss
 * - Auto-recover WebGL context loss
 */

import React, { Suspense, useMemo, useRef, useEffect, useState } from 'react'
import { Canvas, useFrame, useThree } from '@react-three/fiber'
import { OrbitControls } from '@react-three/drei'
import * as THREE from 'three'
import type { Molecule } from '@/lib/molecule'

interface ThreeDViewerProps {
  molecule: Molecule
  selectedAtomId?: string | null
  hoveredAtomId?: string | null
  onAtomHover?: (atomId: string | null) => void
  onAtomClick?: (atomId: string) => void
  width?: number
  height?: number
  className?: string
}

// CPK element colors
const ELEMENT_COLORS: Record<string, number> = {
  H: 0xffffff,
  C: 0x909090,
  N: 0x3050f8,
  O: 0xff0d0d,
  F: 0x90e050,
  P: 0xff8000,
  S: 0xffff30,
  Cl: 0x1ff01f,
  Br: 0xa62929,
  I: 0x9400d3,
  Li: 0xcc80ff,
  Na: 0xab5cf2,
  K: 0x8f40d4,
  Mg: 0x8aff00,
  Ca: 0x3dff00,
  Fe: 0xe06633,
  Cu: 0xc88033,
  Zn: 0x7d80b0,
}

const getElementColor = (element: string): number => {
  return ELEMENT_COLORS[element] || 0xcccccc
}

// Atom sphere component
function AtomSphere({
  atom,
  selected,
  hovered,
  onClick,
  onHover,
}: {
  atom: { id: string; element: string; position: [number, number, number] }
  selected: boolean
  hovered: boolean
  onClick?: () => void
  onHover?: (hovered: boolean) => void
}) {
  const meshRef = useRef<THREE.Mesh>(null)
  const color = getElementColor(atom.element)
  const scale = hovered ? 1.3 : selected ? 1.2 : 1.0
  const emissive = hovered ? 0x444444 : selected ? 0x000044 : 0x000000

  // Animate scale
  useFrame(() => {
    if (meshRef.current) {
      meshRef.current.scale.lerp(new THREE.Vector3(scale, scale, scale), 0.1)
    }
  })

  return (
    <mesh
      ref={meshRef}
      position={atom.position}
      onClick={onClick}
      onPointerOver={() => onHover?.(true)}
      onPointerOut={() => onHover?.(false)}
    >
      <sphereGeometry args={[0.3, 16, 16]} />
      <meshStandardMaterial
        color={color}
        emissive={emissive}
        emissiveIntensity={hovered ? 0.3 : selected ? 0.2 : 0}
        metalness={0.3}
        roughness={0.4}
      />
    </mesh>
  )
}

// Bond cylinder component
function BondCylinder({
  atom1,
  atom2,
  order,
  selected,
}: {
  atom1: { position: [number, number, number] }
  atom2: { position: [number, number, number] }
  order: number
  selected: boolean
}) {
  const [pos1, pos2] = [new THREE.Vector3(...atom1.position), new THREE.Vector3(...atom2.position)]
  const distance = pos1.distanceTo(pos2)
  const midpoint = new THREE.Vector3().addVectors(pos1, pos2).multiplyScalar(0.5)
  const direction = new THREE.Vector3().subVectors(pos2, pos1).normalize()

  // Calculate rotation
  const up = new THREE.Vector3(0, 1, 0)
  const quaternion = new THREE.Quaternion()
  quaternion.setFromUnitVectors(up, direction)

  // Bond radius based on order
  const radius = order === 3 ? 0.08 : order === 2 ? 0.06 : 0.05

  return (
    <mesh position={midpoint} quaternion={quaternion}>
      <cylinderGeometry args={[radius, radius, distance, 16]} />
      <meshStandardMaterial
        color={selected ? 0xffff00 : 0x888888}
        metalness={0.2}
        roughness={0.6}
      />
    </mesh>
  )
}

// Scene component
function MoleculeScene({
  molecule,
  selectedAtomId,
  hoveredAtomId,
  onAtomHover,
  onAtomClick,
}: {
  molecule: Molecule
  selectedAtomId?: string | null
  hoveredAtomId?: string | null
  onAtomHover?: (atomId: string | null) => void
  onAtomClick?: (atomId: string) => void
}) {
  const atoms = molecule.getAtoms()
  const bonds = molecule.getBonds()

  // Create atom lookup
  const atomMap = useMemo(() => {
    const map = new Map<string, { id: string; element: string; position: [number, number, number] }>()
    atoms.forEach(atom => {
      map.set(atom.id, {
        id: atom.id,
        element: atom.element,
        position: atom.position,
      })
    })
    return map
  }, [atoms])

  return (
    <>
      {/* Render bonds */}
      {bonds.map(bond => {
        const atom1 = atomMap.get(bond.atom1)
        const atom2 = atomMap.get(bond.atom2)
        if (!atom1 || !atom2) return null

        return (
          <BondCylinder
            key={bond.id}
            atom1={atom1}
            atom2={atom2}
            order={bond.order}
            selected={selectedAtomId === bond.atom1 || selectedAtomId === bond.atom2}
          />
        )
      })}

      {/* Render atoms */}
      {atoms.map(atom => {
        const atomData = atomMap.get(atom.id)
        if (!atomData) return null

        return (
          <AtomSphere
            key={atom.id}
            atom={atomData}
            selected={selectedAtomId === atom.id}
            hovered={hoveredAtomId === atom.id}
            onClick={() => onAtomClick?.(atom.id)}
            onHover={(hovered) => onAtomHover?.(hovered ? atom.id : null)}
          />
        )
      })}
    </>
  )
}

// WebGL context recovery
function WebGLRecovery() {
  const { gl } = useThree()
  const [lost, setLost] = useState(false)

  useEffect(() => {
    const handleContextLost = (e: Event) => {
      e.preventDefault()
      setLost(true)
      console.warn('WebGL context lost')
    }

    const handleContextRestored = () => {
      setLost(false)
      console.log('WebGL context restored')
    }

    gl.domElement.addEventListener('webglcontextlost', handleContextLost)
    gl.domElement.addEventListener('webglcontextrestored', handleContextRestored)

    return () => {
      gl.domElement.removeEventListener('webglcontextlost', handleContextLost)
      gl.domElement.removeEventListener('webglcontextrestored', handleContextRestored)
    }
  }, [gl])

  if (lost) {
    return (
      <mesh>
        <boxGeometry args={[1, 1, 1]} />
        <meshStandardMaterial color="#ff0000" />
      </mesh>
    )
  }

  return null
}

export function ThreeDViewer({
  molecule,
  selectedAtomId,
  hoveredAtomId,
  onAtomHover,
  onAtomClick,
  width = 400,
  height = 400,
  className = '',
}: ThreeDViewerProps) {
  const [mounted, setMounted] = useState(false)

  // Lazy mount to avoid WebGL loss
  useEffect(() => {
    setMounted(true)
    return () => {
      setMounted(false)
    }
  }, [])

  if (molecule.isEmpty()) {
    return (
      <div
        className={`flex items-center justify-center bg-gray-900 text-gray-500 ${className}`}
        style={{ width, height }}
      >
        <p className="text-sm">No molecule to display</p>
      </div>
    )
  }

  if (!mounted) {
    return (
      <div
        className={`flex items-center justify-center bg-gray-900 text-gray-500 ${className}`}
        style={{ width, height }}
      >
        <p className="text-sm">Loading 3D viewer...</p>
      </div>
    )
  }

  return (
    <ErrorBoundary>
      <div className={className} style={{ width, height }}>
        <Canvas
        camera={{ position: [0, 0, 10], fov: 50 }}
        gl={{ antialias: true, alpha: false, preserveDrawingBuffer: false }}
        style={{ background: '#1a1a1a', width, height }}
      >
        <Suspense fallback={null}>
          <ambientLight intensity={0.6} />
          <directionalLight position={[5, 5, 5]} intensity={0.9} />
          <directionalLight position={[-5, -3, -2]} intensity={0.5} />
          <directionalLight position={[0, 5, -5]} intensity={0.3} />

          <MoleculeScene
            molecule={molecule}
            selectedAtomId={selectedAtomId}
            hoveredAtomId={hoveredAtomId}
            onAtomHover={onAtomHover}
            onAtomClick={onAtomClick}
          />

          <OrbitControls
            enablePan={true}
            enableZoom={true}
            enableRotate={true}
            minDistance={2}
            maxDistance={20}
          />

          <WebGLRecovery />
        </Suspense>
        </Canvas>
      </div>
    </ErrorBoundary>
  )
}

