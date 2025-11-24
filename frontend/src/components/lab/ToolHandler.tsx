import React from 'react'
import { useThree } from '@react-three/fiber'
import { useLabStore } from '../../store/labStore'
import { getTool } from '../../tools'
import { computeAutoBonds } from '../../utils/bondingEngine'
import * as THREE from 'three'

/**
 * ToolHandler routes pointer events to the currently active tool
 */
export default function ToolHandler() {
  const { camera, gl, scene, raycaster } = useThree()
  const currentTool = useLabStore(s => s.currentTool)
  const store = useLabStore.getState()
  const autoBond = useLabStore(s => s.autoBond)

  React.useEffect(() => {
    const handlePointerDown = (event: PointerEvent) => {
      // Only handle if clicking on canvas
      if (event.target !== gl.domElement && !gl.domElement.contains(event.target as Node)) return

      const rect = gl.domElement.getBoundingClientRect()
      const pointer = new THREE.Vector2()
      pointer.x = ((event.clientX - rect.left) / rect.width) * 2 - 1
      pointer.y = -((event.clientY - rect.top) / rect.height) * 2 + 1

      // Raycast to find clicked object
      raycaster.setFromCamera(pointer, camera)
      const intersects = raycaster.intersectObjects(scene.children, true)
      
      const r3fEvent = {
        clientX: event.clientX,
        clientY: event.clientY,
        target: gl.domElement,
        object: intersects.length > 0 ? intersects[0].object : null,
        camera: camera,
      }

      const tool = getTool(currentTool)
      if (tool.onPointerDown) {
        tool.onPointerDown(r3fEvent, store)
        
        // Auto-bond after adding atom if enabled
        if (currentTool === 'add_atom' && autoBond) {
          setTimeout(() => {
            const mol = store.molecule
            const newBonds = computeAutoBonds(mol)
            newBonds.forEach(({ a, b }) => {
              store.addBond(a, b, 1)
            })
          }, 0)
        }
      }
    }

    gl.domElement.addEventListener('pointerdown', handlePointerDown)
    return () => {
      gl.domElement.removeEventListener('pointerdown', handlePointerDown)
    }
  }, [currentTool, autoBond, gl, camera, scene, raycaster, store])

  return null
}

