import React, { Suspense, useEffect, useState } from 'react'
import { Canvas, useThree } from '@react-three/fiber'
import { OrbitControls, ContactShadows, Environment } from '@react-three/drei'
import AtomMesh from './r3f/AtomMesh'
import BondMesh from './r3f/BondMesh'
import { CanvasClickHandler } from './r3f/CanvasClickHandler'
import { useMoleculeStore } from '../store/moleculeStore'
import { moleculeToRenderable, removeAtom, removeBond, updateAtomPosition } from '../lib/engineAdapter'
import { selectionManager } from './r3f/SelectionManager'
import { useBondTool } from './r3f/BondTool'
import { screenToWorld } from '../lib/raycasting'
import { useAutoBondOnDrop } from '../hooks/useAutoBondOnDrop'

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

      // Convert screen to world coordinates using plane at atom's Y position
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

      // Clamp position to reasonable bounds
      const clampedPos: [number, number, number] = [
        Math.max(-5, Math.min(5, worldPos.x)),
        Math.max(-5, Math.min(5, worldPos.y)),
        Math.max(-5, Math.min(5, worldPos.z)),
      ]

      // Update atom position
      updateAtomPosition(dragAtomId, clampedPos)
    }

    window.addEventListener('pointermove', handlePointerMove)
    return () => {
      window.removeEventListener('pointermove', handlePointerMove)
    }
  }, [isDragging, dragAtomId, currentMolecule, camera, gl])

  // Invisible plane for dragging (not needed but kept for reference)
  return null
}

// Canvas click handler - handled in main component

export default function MoleculeViewer() {
  const currentMolecule = useMoleculeStore((state) => state.currentMolecule)
  const fetchPredictions = useMoleculeStore((state) => state.fetchPredictions)
  const tool = useMoleculeStore((state) => state.tool)
  const selectedAtomId = useMoleculeStore((state) => state.selectedAtomId)
  const selectedBondId = useMoleculeStore((state) => state.selectedBondId)

  const { atoms, bonds } = moleculeToRenderable(currentMolecule)

  // Use bond tool
  useBondTool()
  
  // Auto-bond on atom drop (when autoBond is enabled)
  useAutoBondOnDrop()

  // Auto-update predictions when molecule changes
  useEffect(() => {
    if (currentMolecule) {
      fetchPredictions()
    }
  }, [currentMolecule, fetchPredictions])

  // Track cursor state
  const [cursor, setCursor] = useState<'default' | 'pointer' | 'grabbing' | 'crosshair'>('default')

  useEffect(() => {
    const updateCursor = () => {
      const draggingId = selectionManager.getDraggingAtomId()
      const hoveredId = selectionManager.getHoveredAtomId()
      
      if (draggingId) {
        setCursor('grabbing')
      } else if (tool === 'add-atom') {
        setCursor('crosshair')
      } else if (hoveredId) {
        setCursor('pointer')
      } else {
        setCursor('default')
      }
    }

    const unsubHover = selectionManager.on('hover', updateCursor)
    const unsubDrag = selectionManager.on('drag', updateCursor)

    updateCursor() // Initial update

    return () => {
      unsubHover()
      unsubDrag()
    }
  }, [tool])

  // Handle delete key
  useEffect(() => {
    const handleKeyDown = (e: KeyboardEvent) => {
      if ((e.key === 'Delete' || e.key === 'Backspace') && tool === 'delete') {
        if (selectedAtomId) {
          removeAtom(selectedAtomId)
        } else if (selectedBondId) {
          removeBond(selectedBondId)
        }
      }
    }

    window.addEventListener('keydown', handleKeyDown)
    return () => window.removeEventListener('keydown', handleKeyDown)
  }, [tool, selectedAtomId, selectedBondId])

  return (
    <div className="relative w-full h-full">
      <div
        style={{ cursor }}
        className="w-full h-full rounded-lg overflow-hidden bg-spaceGrey"
      >
        <Canvas dpr={[1, 2]} camera={{ position: [0, 0, 12], fov: 45 }} style={{ width: '100%', height: '100%' }}>
          {/* Cold white/ivory lighting - NO blue */}
          <ambientLight intensity={0.4} color="#F6F7F8" />
          <directionalLight position={[5, 10, 7]} intensity={1.2} color="#FFFFFF" castShadow />
          <hemisphereLight intensity={0.35} groundColor={0x2B2E33 as any} skyColor={0xF6F7F8 as any} />
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
              <BondMesh 
                key={bond.id} 
                id={bond.id}
                from={bond.from} 
                to={bond.to}
                order={bond.order}
              />
            ))}
            <InteractionLayer />
            <CanvasClickHandler />
            <ContactShadows opacity={0.2} blur={2.5} far={12} resolution={256} color="#0F1115" />
        </Suspense>
          {/* Chrome-like environment with cold white reflections */}
          <Environment preset="city" environmentIntensity={0.8} />
          <OrbitControls
            enablePan={tool === 'select'} 
            enableZoom={true} 
            enableRotate={true}
            enableDamping
            dampingFactor={0.08}
          />
      </Canvas>
      </div>
      
      {/* Property badges */}
      {(selectedAtomId || selectedBondId) && (
        <PropertyBadge 
          atomId={selectedAtomId} 
          bondId={selectedBondId}
          molecule={currentMolecule}
        />
      )}
    </div>
  )
}

// Property badge component
function PropertyBadge({ 
  atomId, 
  bondId, 
  molecule 
}: { 
  atomId: string | null
  bondId: string | null
  molecule: any
}) {
  if (!molecule) return null

  if (atomId) {
    const atom = molecule.atoms.get(atomId)
    if (!atom) return null

    return (
      <div className="absolute bottom-4 left-4 frosted-glass rounded-lg shadow-glass border border-chrome/20 p-3 text-sm backdrop-blur-md">
        <div className="font-semibold text-ivory mb-1">Atom: {atom.element}</div>
        <div className="text-chrome">
          Position: ({atom.position[0].toFixed(2)}, {atom.position[1].toFixed(2)}, {atom.position[2].toFixed(2)})
        </div>
        <div className="text-chrome/70">
          ID: {atom.id}
        </div>
      </div>
    )
  }

  if (bondId) {
    const bond = molecule.bonds.get(bondId)
    if (!bond) return null

    const atom1 = molecule.atoms.get(bond.a1)
    const atom2 = molecule.atoms.get(bond.a2)
    if (!atom1 || !atom2) return null

    const dx = atom2.position[0] - atom1.position[0]
    const dy = atom2.position[1] - atom1.position[1]
    const dz = atom2.position[2] - atom1.position[2]
    const length = Math.sqrt(dx * dx + dy * dy + dz * dz)

    return (
      <div className="absolute bottom-4 left-4 frosted-glass rounded-lg shadow-glass border border-chrome/20 p-3 text-sm backdrop-blur-md">
        <div className="font-semibold text-ivory mb-1">
          Bond: {bond.order === 1 ? 'Single' : bond.order === 2 ? 'Double' : 'Triple'}
        </div>
        <div className="text-chrome">
          Length: {length.toFixed(2)} Å
        </div>
        <div className="text-chrome">
          Between: {atom1.element} - {atom2.element}
        </div>
      </div>
    )
  }

  return null
}

