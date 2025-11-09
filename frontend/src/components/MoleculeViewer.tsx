import React, { Suspense, useEffect, useState } from 'react'
import { Canvas, useThree } from '@react-three/fiber'
import { OrbitControls } from '@react-three/drei'
import AtomMesh from './r3f/AtomMesh'
import BondMesh from './r3f/BondMesh'
import { useMoleculeStore } from '../store/moleculeStore'
import { moleculeToRenderable, updateAtomPosition } from '../lib/engineAdapter'
import { selectionManager } from './r3f/SelectionManager'
import { useBondTool } from './r3f/BondTool'
import { screenToWorld } from '../lib/raycasting'

// Interaction layer component
function InteractionLayer() {
  const { camera, gl } = useThree()
  const currentMolecule = useMoleculeStore((state) => state.currentMolecule)
  const [isDragging, setIsDragging] = useState(false)
  const [dragAtomId, setDragAtomId] = useState<string | null>(null)

  // Handle drag updates
  useEffect(() => {
    const handleDrag = () => {
      const draggingId = selectionManager.getDraggingAtomId()
      setIsDragging(draggingId !== null)
      setDragAtomId(draggingId)
    }

    const unsubscribe = selectionManager.on('drag', handleDrag)
    return unsubscribe
  }, [])

  // Handle pointer move for dragging
  useEffect(() => {
    if (!isDragging || !dragAtomId || !currentMolecule) return

    const handlePointerMove = (e: PointerEvent) => {
      const atom = currentMolecule.atoms.get(dragAtomId)
      if (!atom) return

      // Get canvas element for coordinate conversion
      const canvas = gl.domElement
      if (!canvas) return

      // Convert screen to world coordinates
      const worldPos = screenToWorld(
        {
          ...e,
          target: canvas,
          clientX: e.clientX,
          clientY: e.clientY,
        } as any,
        camera,
        atom.position[1] // Use atom's Y position as plane
      )

      // Update atom position
      updateAtomPosition(dragAtomId, [worldPos.x, worldPos.y, worldPos.z])
    }

    window.addEventListener('pointermove', handlePointerMove)
    return () => {
      window.removeEventListener('pointermove', handlePointerMove)
    }
  }, [isDragging, dragAtomId, currentMolecule, camera, gl])

  // Invisible plane for dragging (not needed but kept for reference)
  return null
}

export default function MoleculeViewer() {
  const currentMolecule = useMoleculeStore((state) => state.currentMolecule)
  const fetchPredictions = useMoleculeStore((state) => state.fetchPredictions)

  const { atoms, bonds } = moleculeToRenderable(currentMolecule)

  // Use bond tool
  useBondTool()

  // Auto-update predictions when molecule changes
  useEffect(() => {
    if (currentMolecule) {
      fetchPredictions()
    }
  }, [currentMolecule, fetchPredictions])

  // Track cursor state
  const [cursor, setCursor] = useState<'default' | 'pointer' | 'grabbing'>('default')

  useEffect(() => {
    const updateCursor = () => {
      const draggingId = selectionManager.getDraggingAtomId()
      const hoveredId = selectionManager.getHoveredAtomId()
      
      if (draggingId) {
        setCursor('grabbing')
      } else if (hoveredId) {
        setCursor('pointer')
      } else {
        setCursor('default')
      }
    }

    const unsubHover = selectionManager.on('hover', updateCursor)
    const unsubDrag = selectionManager.on('drag', updateCursor)

    return () => {
      unsubHover()
      unsubDrag()
    }
  }, [])

  return (
    <div
      style={{ height: 420, cursor }}
      className="rounded-lg overflow-hidden bg-aluminum-light"
    >
      <Canvas dpr={[1, 2]} camera={{ position: [0, 0, 12], fov: 45 }}>
        <ambientLight intensity={0.6} />
        <directionalLight position={[5, 10, 7]} intensity={0.6} />
        <Suspense fallback={null}>
          {atoms.map((atom) => (
            <AtomMesh
              key={atom.id}
              id={atom.id}
              position={atom.position}
              element={atom.element as any}
            />
          ))}
          {bonds.map((bond) => (
            <BondMesh key={bond.id} from={bond.from} to={bond.to} />
          ))}
          <InteractionLayer />
        </Suspense>
        <OrbitControls enablePan={true} enableZoom={true} enableRotate={true} />
      </Canvas>
    </div>
  )
}

