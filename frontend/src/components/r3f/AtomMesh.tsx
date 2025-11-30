import React, { useRef, useState } from 'react'
import * as THREE from 'three'
import { useFrame } from '@react-three/fiber'
import { Outlines, Text } from '@react-three/drei'
import { selectionManager } from './SelectionManager'
import { useMoleculeStore } from '../../store/moleculeStore'
import { kernelSelectionManager } from '../../kernel/selectionManager'
import type { ColorScheme, RenderMode } from '../../store/moleculeStore'
import { ELEMENT_ACCENT_COLORS, ELEMENT_BASE_COLORS, ELEMENT_RADII } from '../../utils/chemistry'

interface AtomMeshProps {
  id: string
  position: [number, number, number]
  element: 'C' | 'H' | 'O' | 'N' | 'F' | 'S' | 'P' | 'Cl' | 'Br' | 'I'
  renderMode: RenderMode
  colorScheme: ColorScheme
  colorOverride?: number
  highlighted?: boolean
}

export default function AtomMesh({
  id,
  position,
  element,
  renderMode,
  colorScheme,
  colorOverride,
  highlighted = false,
}: AtomMeshProps) {
  const meshRef = useRef<THREE.Mesh>(null)
  const [isHovered, setIsHovered] = useState(false)
  const [isDragging, setIsDragging] = useState(false)
  const selectedAtomId = useMoleculeStore((state) => state.selectedAtomId)
  const tool = useMoleculeStore((state) => state.tool)
  
  const radius = ELEMENT_RADII[element] || 1.0
  const isHighlighted = highlighted
  const baseColorBase = colorOverride ?? ELEMENT_BASE_COLORS[element] ?? 0x9da3ae
  const baseColor = isHighlighted ? 0xffc857 : baseColorBase
  const accentColor = ELEMENT_ACCENT_COLORS[element] || 0x2B2E33
  const isSelected = selectedAtomId === id
  
  // Visual feedback based on tool - using neonCyan for selection
  let scale = isSelected ? 1.15 : 1.0
  if (renderMode === 'spacefill') {
    scale *= 1.35
  } else if (renderMode === 'wireframe') {
    scale *= 0.75
  }
  if (isHighlighted) {
    scale *= 1.1
  }
  let outlineColor = 0x8BF3FF // neonCyan for select/hover
  if (tool === 'delete' && (isHovered || isSelected)) {
    outlineColor = 0xC6BDFE // violetEdge for delete
  } else if (tool === 'bond' && (isHovered || isSelected)) {
    outlineColor = 0x8BF3FF // neonCyan for bond
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
      // Update both UI and kernel selection managers
      selectionManager.onSelect(id)
      kernelSelectionManager.selectAtom(id)
      useMoleculeStore.getState().selectAtom(id)
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

  const segments = (typeof window !== 'undefined' && window.devicePixelRatio && window.devicePixelRatio > 1.5) ? 48 : 64

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
      <sphereGeometry args={[radius, segments, segments]} />
      {renderMode === 'wireframe' ? (
        <meshBasicMaterial color={baseColor} wireframe transparent opacity={0.9} />
      ) : (
        <meshPhysicalMaterial
          color={baseColor}
          roughness={renderMode === 'spacefill' ? 0.35 : 0.2}
          metalness={renderMode === 'spacefill' ? 0.2 : 0.4}
          clearcoat={1}
          clearcoatRoughness={0.1}
          thickness={0.5}
          transparent={isDragging}
          opacity={opacity}
          envMapIntensity={1.2}
        />
      )}
      {renderMode !== 'wireframe' && (isHovered || isSelected || isHighlighted) && (
        <>
          <Outlines thickness={0.08} color={accentColor} />
          <Outlines thickness={0.12} color={outlineColor} />
        </>
      )}
      {renderMode === 'ballstick' && (
        <Text
          position={[0, 0, radius * 1.05]}
          fontSize={radius * 0.5}
          color={isSelected || isHovered ? "#000000" : "#2B2E33"}
          anchorX="center"
          anchorY="middle"
          outlineWidth={0.01}
          outlineColor="#FFFFFF"
          renderOrder={1000}
        >
          {element}
        </Text>
      )}
    </mesh>
  )
}

