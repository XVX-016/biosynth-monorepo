/**
 * Editor2D - 2D molecule editor using Konva
 * 
 * Handles atom placement, bond drawing, selection, and grid snapping
 */

import React, { useRef, useEffect, useState, useCallback } from 'react'
import { Stage, Layer, Circle, Line, Text } from 'react-konva'
import { moleculeEngine } from '../engines/MoleculeStateEngine'
import { toolController, ToolMode, type Point } from '../engines/ToolController'
import { commandStack, AddAtomCommand, AddBondCommand, RemoveAtomCommand, MoveAtomCommand } from '../engines/CommandStack'

interface Editor2DProps {
  width?: number
  height?: number
  gridSize?: number
  selectedElement?: string | null
  onAtomSelect?: (atomId: string | null) => void
}

const HITBOX_RADIUS = 20 // Pixels
const GRID_SIZE = 40 // Default grid size

// Element colors (CPK)
const ELEMENT_COLORS: Record<string, string> = {
  H: '#FFFFFF',
  C: '#909090',
  N: '#3050F8',
  O: '#FF0D0D',
  F: '#90E050',
  Cl: '#1FF01F',
  Br: '#A62929',
  I: '#940094',
  S: '#FFFF30',
  P: '#FF8000',
}

export default function Editor2D({
  width = 800,
  height = 600,
  gridSize = GRID_SIZE,
  selectedElement = null,
  onAtomSelect,
}: Editor2DProps) {
  const stageRef = useRef<any>(null)
  const [draggingAtomId, setDraggingAtomId] = useState<string | null>(null)
  const [bondPreview, setBondPreview] = useState<{ start: Point; end: Point } | null>(null)

  // Snap point to grid
  const snapToGrid = useCallback((x: number, y: number): Point => {
    return {
      x: Math.round(x / gridSize) * gridSize,
      y: Math.round(y / gridSize) * gridSize,
    }
  }, [gridSize])

  // Find atom at point (hitbox detection)
  const findAtomAtPoint = useCallback((point: Point): string | null => {
    const atoms = moleculeEngine.getAllAtoms()
    let closest: { id: string; dist: number } | null = null

    atoms.forEach(atom => {
      const dx = atom.x - point.x
      const dy = atom.y - point.y
      const dist = Math.sqrt(dx * dx + dy * dy)
      
      if (dist < HITBOX_RADIUS) {
        if (!closest || dist < closest.dist) {
          closest = { id: atom.id, dist }
        }
      }
    })

    return closest?.id || null
  }, [])

  // Handle mouse down
  const handleMouseDown = useCallback((e: any) => {
    const stage = e.target.getStage()
    const pointerPos = stage.getPointerPosition()
    const snappedPos = snapToGrid(pointerPos.x, pointerPos.y)
    const hitAtomId = findAtomAtPoint(snappedPos)

    if (toolController.activeTool === ToolMode.addAtom && selectedElement) {
      // Add atom
      const cmd = new AddAtomCommand(selectedElement, snappedPos.x, snappedPos.y)
      commandStack.run(cmd)
      onAtomSelect?.(null)
    } else if (toolController.activeTool === ToolMode.addBond) {
      // Bond drawing logic
      if (hitAtomId) {
        if (toolController.bondPreviewStart === null) {
          // First click - start bond
          toolController.bondPreviewStart = hitAtomId
        } else if (toolController.bondPreviewStart !== hitAtomId) {
          // Second click - complete bond
          const cmd = new AddBondCommand(toolController.bondPreviewStart, hitAtomId, 1)
          commandStack.run(cmd)
          toolController.bondPreviewStart = null
          setBondPreview(null)
        }
      }
    } else if (toolController.activeTool === ToolMode.erase && hitAtomId) {
      // Erase atom
      const cmd = new RemoveAtomCommand(hitAtomId)
      commandStack.run(cmd)
    } else if (toolController.activeTool === ToolMode.select || toolController.activeTool === ToolMode.cursor) {
      // Select atom
      toolController.handleMouseDown(snappedPos, hitAtomId)
      onAtomSelect?.(hitAtomId)
      if (hitAtomId) {
        setDraggingAtomId(hitAtomId)
      }
    }
  }, [selectedElement, snapToGrid, findAtomAtPoint, onAtomSelect])

  // Handle mouse move
  const handleMouseMove = useCallback((e: any) => {
    const stage = e.target.getStage()
    const pointerPos = stage.getPointerPosition()

    // Update bond preview
    if (toolController.activeTool === ToolMode.addBond && toolController.bondPreviewStart) {
      const startAtom = moleculeEngine.getAtom(toolController.bondPreviewStart)
      if (startAtom) {
        setBondPreview({
          start: { x: startAtom.x, y: startAtom.y },
          end: pointerPos,
        })
      }
    } else {
      setBondPreview(null)
    }

    // Handle dragging
    if (draggingAtomId) {
      const snappedPos = snapToGrid(pointerPos.x, pointerPos.y)
      const atom = moleculeEngine.getAtom(draggingAtomId)
      if (atom) {
        const cmd = new MoveAtomCommand(draggingAtomId, snappedPos.x, snappedPos.y)
        commandStack.run(cmd)
      }
    }
  }, [draggingAtomId, snapToGrid])

  // Handle mouse up
  const handleMouseUp = useCallback(() => {
    setDraggingAtomId(null)
    toolController.handleMouseUp({ x: 0, y: 0 }, null)
  }, [])

  // Render atoms
  const renderAtoms = () => {
    return moleculeEngine.getAllAtoms().map(atom => {
      const isSelected = toolController.selectedAtom === atom.id
      const color = ELEMENT_COLORS[atom.element] || '#909090'
      
      return (
        <React.Fragment key={atom.id}>
          <Circle
            x={atom.x}
            y={atom.y}
            radius={16}
            fill={color}
            stroke={isSelected ? '#FFD700' : '#000'}
            strokeWidth={isSelected ? 3 : 1}
            onClick={() => {
              toolController.selectedAtom = atom.id
              onAtomSelect?.(atom.id)
            }}
          />
          <Text
            x={atom.x - 6}
            y={atom.y - 8}
            text={atom.element}
            fontSize={12}
            fill="#000"
            fontStyle="bold"
          />
        </React.Fragment>
      )
    })
  }

  // Render bonds
  const renderBonds = () => {
    return moleculeEngine.getAllBonds().map(bond => {
      const atomA = moleculeEngine.getAtom(bond.atoms[0])
      const atomB = moleculeEngine.getAtom(bond.atoms[1])
      
      if (!atomA || !atomB) return null

      return (
        <Line
          key={bond.id}
          points={[atomA.x, atomA.y, atomB.x, atomB.y]}
          stroke="#000"
          strokeWidth={bond.order === 1 ? 2 : bond.order === 2 ? 3 : 4}
        />
      )
    })
  }

  // Render grid
  const renderGrid = () => {
    const lines = []
    for (let i = 0; i <= width; i += gridSize) {
      lines.push(
        <Line
          key={`v-${i}`}
          points={[i, 0, i, height]}
          stroke="#E0E0E0"
          strokeWidth={0.5}
        />
      )
    }
    for (let i = 0; i <= height; i += gridSize) {
      lines.push(
        <Line
          key={`h-${i}`}
          points={[0, i, width, i]}
          stroke="#E0E0E0"
          strokeWidth={0.5}
        />
      )
    }
    return lines
  }

  return (
    <div className="w-full h-full bg-white border border-gray-300">
      <Stage
        ref={stageRef}
        width={width}
        height={height}
        onMouseDown={handleMouseDown}
        onMouseMove={handleMouseMove}
        onMouseUp={handleMouseUp}
      >
        <Layer>
          {/* Grid */}
          {renderGrid()}
          
          {/* Bonds */}
          {renderBonds()}
          
          {/* Bond preview */}
          {bondPreview && (
            <Line
              points={[bondPreview.start.x, bondPreview.start.y, bondPreview.end.x, bondPreview.end.y]}
              stroke="#888"
              strokeWidth={2}
              dash={[5, 5]}
            />
          )}
          
          {/* Atoms */}
          {renderAtoms()}
        </Layer>
      </Stage>
    </div>
  )
}

