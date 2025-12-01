/**
 * CanvasLayer - Pure 2D canvas rendering for molecule editor
 * 
 * Phase 2: Rewrite Drawing Layer
 * 
 * Features:
 * - Atoms first, bonds second (proper z-ordering)
 * - Hover highlights
 * - Selection outlines
 * - Proper hitboxes (8px radius)
 * - Pixel ratio fix
 * - No coordinate drift
 */

import React, { useRef, useEffect, useCallback, useMemo } from 'react'
import type { Molecule } from '@/lib/molecule'
import type { AtomImpl, BondImpl } from '@/lib/molecule'
import { ELEMENT_COLORS } from './constants'

import type { PointerManager } from '@/lib/molecule/input'

interface CanvasLayerProps {
  molecule: Molecule | null
  selectedAtomId: string | null
  selectedBondId: string | null
  hoveredAtomId: string | null
  hoveredBondId: string | null
  width: number
  height: number
  onAtomClick?: (atomId: string, event: MouseEvent) => void
  onBondClick?: (bondId: string, event: MouseEvent) => void
  onAtomHover?: (atomId: string | null) => void
  onBondHover?: (bondId: string | null) => void
  scale?: number
  offsetX?: number
  offsetY?: number
  pointerManager?: PointerManager
  bondStartAtomId?: string | null
}

// Element colors for rendering
const ELEMENT_COLORS: Record<string, string> = {
  H: '#ffffff',
  C: '#909090',
  N: '#3050f8',
  O: '#ff0d0d',
  S: '#ffff30',
  P: '#ff8000',
  F: '#90e050',
  Cl: '#1ff01f',
  Br: '#a62929',
  I: '#940094',
}

// Hit detection radius (8px as specified)
const HIT_RADIUS = 8

export function CanvasLayer({
  molecule,
  selectedAtomId,
  selectedBondId,
  hoveredAtomId,
  hoveredBondId,
  width,
  height,
  onAtomClick,
  onBondClick,
  onAtomHover,
  onBondHover,
  scale = 1,
  offsetX = 0,
  offsetY = 0,
  pointerManager,
  bondStartAtomId,
}: CanvasLayerProps) {
  const canvasRef = useRef<HTMLCanvasElement>(null)
  const rafRef = useRef<number>()

  // Convert 3D to 2D projection (orthographic, looking down Z-axis)
  const projectTo2D = useCallback((position: [number, number, number]): [number, number] => {
    const [x, y] = position
    return [
      (x + offsetX) * scale + width / 2,
      (y + offsetY) * scale + height / 2,
    ]
  }, [width, height, scale, offsetX, offsetY])

  // Draw atom
  const drawAtom = useCallback((
    ctx: CanvasRenderingContext2D,
    atom: AtomImpl,
    isSelected: boolean,
    isHovered: boolean,
  ) => {
    const [x, y] = projectTo2D(atom.position)
    const element = atom.element
    const color = ELEMENT_COLORS[element] || '#909090'

    // Atom circle
    ctx.beginPath()
    ctx.arc(x, y, 12, 0, Math.PI * 2)
    ctx.fillStyle = color
    ctx.fill()

    // Border
    ctx.strokeStyle = isSelected ? '#3b82f6' : isHovered ? '#60a5fa' : '#000000'
    ctx.lineWidth = isSelected ? 3 : isHovered ? 2 : 1
    ctx.stroke()

    // Element label
    ctx.fillStyle = '#000000'
    ctx.font = '12px Arial'
    ctx.textAlign = 'center'
    ctx.textBaseline = 'middle'
    ctx.fillText(element, x, y)

    // Selection outline
    if (isSelected) {
      ctx.beginPath()
      ctx.arc(x, y, 18, 0, Math.PI * 2)
      ctx.strokeStyle = '#3b82f6'
      ctx.lineWidth = 2
      ctx.setLineDash([5, 5])
      ctx.stroke()
      ctx.setLineDash([])
    }
  }, [projectTo2D])

  // Draw bond
  const drawBond = useCallback((
    ctx: CanvasRenderingContext2D,
    bond: BondImpl,
    atom1: AtomImpl,
    atom2: AtomImpl,
    isSelected: boolean,
    isHovered: boolean,
  ) => {
    const [x1, y1] = projectTo2D(atom1.position)
    const [x2, y2] = projectTo2D(atom2.position)

    // Bond line
    ctx.beginPath()
    ctx.moveTo(x1, y1)
    ctx.lineTo(x2, y2)
    ctx.strokeStyle = isSelected ? '#3b82f6' : isHovered ? '#60a5fa' : '#666666'
    ctx.lineWidth = isSelected ? 3 : isHovered ? 2.5 : bond.order === 1 ? 2 : bond.order === 2 ? 3 : 4
    ctx.stroke()

    // Double/triple bond lines
    if (bond.order === 2) {
      const dx = x2 - x1
      const dy = y2 - y1
      const len = Math.sqrt(dx * dx + dy * dy)
      const perpX = -dy / len * 3
      const perpY = dx / len * 3

      ctx.beginPath()
      ctx.moveTo(x1 + perpX, y1 + perpY)
      ctx.lineTo(x2 + perpX, y2 + perpY)
      ctx.stroke()
    } else if (bond.order === 3) {
      const dx = x2 - x1
      const dy = y2 - y1
      const len = Math.sqrt(dx * dx + dy * dy)
      const perpX = -dy / len * 4
      const perpY = dx / len * 4

      ctx.beginPath()
      ctx.moveTo(x1 + perpX, y1 + perpY)
      ctx.lineTo(x2 + perpX, y2 + perpY)
      ctx.stroke()

      ctx.beginPath()
      ctx.moveTo(x1 - perpX, y1 - perpY)
      ctx.lineTo(x2 - perpX, y2 - perpY)
      ctx.stroke()
    }
  }, [projectTo2D])

  // Render function
  const render = useCallback(() => {
    const canvas = canvasRef.current
    if (!canvas) return

    const ctx = canvas.getContext('2d')
    if (!ctx) return

    // Clear canvas
    ctx.clearRect(0, 0, width, height)

    if (!molecule || molecule.isEmpty()) {
      return
    }

    // Get all atoms and bonds
    const atoms = molecule.getAtoms()
    const bonds = molecule.getBonds()

    // Render bonds first (so atoms appear on top)
    bonds.forEach(bond => {
      const atom1 = molecule.getAtom(bond.atom1)
      const atom2 = molecule.getAtom(bond.atom2)
      if (!atom1 || !atom2) return

      const isSelected = bond.id === selectedBondId
      const isHovered = bond.id === hoveredBondId
      drawBond(ctx, bond, atom1, atom2, isSelected, isHovered)
    })

    // Render atoms second (on top of bonds)
    atoms.forEach(atom => {
      const isSelected = atom.id === selectedAtomId
      const isHovered = atom.id === hoveredAtomId
      drawAtom(ctx, atom, isSelected, isHovered)
    })

    // Draw bond creation preview (if in bond mode with start atom)
    if (bondStartAtomId && hoveredAtomId && hoveredAtomId !== bondStartAtomId) {
      const startAtom = molecule.getAtom(bondStartAtomId)
      const endAtom = molecule.getAtom(hoveredAtomId)
      if (startAtom && endAtom) {
        const [x1, y1] = projectTo2D(startAtom.position)
        const [x2, y2] = projectTo2D(endAtom.position)
        ctx.beginPath()
        ctx.moveTo(x1, y1)
        ctx.lineTo(x2, y2)
        ctx.strokeStyle = '#60a5fa'
        ctx.lineWidth = 2
        ctx.setLineDash([5, 5])
        ctx.stroke()
        ctx.setLineDash([])
      }
    }
  }, [molecule, selectedAtomId, selectedBondId, hoveredAtomId, hoveredBondId, bondStartAtomId, width, height, drawAtom, drawBond, projectTo2D])

  // Handle pixel ratio
  useEffect(() => {
    const canvas = canvasRef.current
    if (!canvas) return

    const dpr = window.devicePixelRatio || 1
    const rect = canvas.getBoundingClientRect()

    // Set actual size in memory (scaled for DPR)
    canvas.width = rect.width * dpr
    canvas.height = rect.height * dpr

    // Scale context to match DPR
    const ctx = canvas.getContext('2d')
    if (ctx) {
      ctx.scale(dpr, dpr)
    }

    // Set display size (CSS pixels)
    canvas.style.width = `${rect.width}px`
    canvas.style.height = `${rect.height}px`

    render()
  }, [width, height, render])

  // Re-render when molecule or selection changes
  useEffect(() => {
    if (rafRef.current) {
      cancelAnimationFrame(rafRef.current)
    }
    rafRef.current = requestAnimationFrame(render)
    return () => {
      if (rafRef.current) {
        cancelAnimationFrame(rafRef.current)
      }
    }
  }, [render])

  // Hit detection: find atom at point
  const findAtomAtPoint = useCallback((x: number, y: number): string | null => {
    if (!molecule) return null

    const atoms = molecule.getAtoms()
    for (const atom of atoms) {
      const [ax, ay] = projectTo2D(atom.position)
      const dx = x - ax
      const dy = y - ay
      const dist = Math.sqrt(dx * dx + dy * dy)
      if (dist <= HIT_RADIUS) {
        return atom.id
      }
    }
    return null
  }, [molecule, projectTo2D])

  // Hit detection: find bond at point
  const findBondAtPoint = useCallback((x: number, y: number): string | null => {
    if (!molecule) return null

    const bonds = molecule.getBonds()
    for (const bond of bonds) {
      const atom1 = molecule.getAtom(bond.atom1)
      const atom2 = molecule.getAtom(bond.atom2)
      if (!atom1 || !atom2) continue

      const [x1, y1] = projectTo2D(atom1.position)
      const [x2, y2] = projectTo2D(atom2.position)

      // Distance from point to line segment
      const dx = x2 - x1
      const dy = y2 - y1
      const len = Math.sqrt(dx * dx + dy * dy)
      if (len === 0) continue

      const t = Math.max(0, Math.min(1, ((x - x1) * dx + (y - y1) * dy) / (len * len)))
      const projX = x1 + t * dx
      const projY = y1 + t * dy
      const dist = Math.sqrt((x - projX) ** 2 + (y - projY) ** 2)

      if (dist <= HIT_RADIUS) {
        return bond.id
      }
    }
    return null
  }, [molecule, projectTo2D])

  // Setup pointer manager integration
  useEffect(() => {
    if (!pointerManager) return

    const canvas = canvasRef.current
    if (!canvas) return

    const handlePointerDown = (e: MouseEvent | TouchEvent) => {
      pointerManager.processEvent(e, 'down')
    }

    const handlePointerMove = (e: MouseEvent | TouchEvent) => {
      pointerManager.processEvent(e, 'move')
    }

    const handlePointerUp = (e: MouseEvent | TouchEvent) => {
      pointerManager.processEvent(e, 'up')
    }

    canvas.addEventListener('mousedown', handlePointerDown)
    canvas.addEventListener('mousemove', handlePointerMove)
    canvas.addEventListener('mouseup', handlePointerUp)
    canvas.addEventListener('touchstart', handlePointerDown)
    canvas.addEventListener('touchmove', handlePointerMove)
    canvas.addEventListener('touchend', handlePointerUp)

    return () => {
      canvas.removeEventListener('mousedown', handlePointerDown)
      canvas.removeEventListener('mousemove', handlePointerMove)
      canvas.removeEventListener('mouseup', handlePointerUp)
      canvas.removeEventListener('touchstart', handlePointerDown)
      canvas.removeEventListener('touchmove', handlePointerMove)
      canvas.removeEventListener('touchend', handlePointerUp)
    }
  }, [pointerManager])

  // Mouse event handlers (for compatibility)
  const handleMouseMove = useCallback((e: React.MouseEvent<HTMLCanvasElement>) => {
    const canvas = canvasRef.current
    if (!canvas) return

    const rect = canvas.getBoundingClientRect()
    const x = e.clientX - rect.left
    const y = e.clientY - rect.top

    // Check atoms first (they're on top)
    const atomId = findAtomAtPoint(x, y)
    if (atomId) {
      onAtomHover?.(atomId)
      onBondHover?.(null)
      return
    }

    // Then check bonds
    const bondId = findBondAtPoint(x, y)
    if (bondId) {
      onBondHover?.(bondId)
      onAtomHover?.(null)
      return
    }

    // Nothing hovered
    onAtomHover?.(null)
    onBondHover?.(null)
  }, [findAtomAtPoint, findBondAtPoint, onAtomHover, onBondHover])

  const handleMouseClick = useCallback((e: React.MouseEvent<HTMLCanvasElement>) => {
    const canvas = canvasRef.current
    if (!canvas) return

    const rect = canvas.getBoundingClientRect()
    const x = e.clientX - rect.left
    const y = e.clientY - rect.top

    // Check atoms first
    const atomId = findAtomAtPoint(x, y)
    if (atomId) {
      onAtomClick?.(atomId, e.nativeEvent)
      return
    }

    // Then check bonds
    const bondId = findBondAtPoint(x, y)
    if (bondId) {
      onBondClick?.(bondId, e.nativeEvent)
      return
    }
  }, [findAtomAtPoint, findBondAtPoint, onAtomClick, onBondClick])

  const handleMouseLeave = useCallback(() => {
    onAtomHover?.(null)
    onBondHover?.(null)
    pointerManager?.cancel()
  }, [onAtomHover, onBondHover, pointerManager])

  return (
    <canvas
      ref={canvasRef}
      width={width}
      height={height}
      onMouseMove={handleMouseMove}
      onClick={handleMouseClick}
      onMouseLeave={handleMouseLeave}
      style={{
        display: 'block',
        cursor: hoveredAtomId || hoveredBondId ? 'pointer' : 'default',
      }}
    />
  )
}

