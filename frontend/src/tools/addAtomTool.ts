import type { Tool } from './toolInterface'
import { screenToWorld } from '../lib/raycasting'

/**
 * addAtomTool expects pointer events from the canvas.
 * Uses screenToWorld to convert click position to 3D coordinates on y=0 plane.
 */
const addAtomTool: Tool = {
  name: 'add_atom',
  onPointerDown: (ev: any, store: any) => {
    // ev: three pointer event from r3f; camera available at ev.camera
    if (!ev.camera || !ev.target) return
    
    // Use existing screenToWorld utility
    const worldPos = screenToWorld(
      {
        clientX: ev.clientX,
        clientY: ev.clientY,
        target: ev.target,
      } as any,
      ev.camera,
      0 // y=0 plane
    )
    
    const pos: [number, number, number] = [worldPos.x, worldPos.y, worldPos.z]
    store.addAtom(store.currentElement, pos)
  }
}

export default addAtomTool

