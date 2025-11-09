import * as THREE from 'three'
import { Raycaster, Vector3, Intersection } from 'three'

/**
 * Get atom under cursor using raycasting
 */
export function getAtomUnderCursor(
  event: MouseEvent | React.PointerEvent,
  camera: THREE.Camera,
  scene: THREE.Scene,
  atomMeshes: THREE.Mesh[]
): THREE.Mesh | null {
  const raycaster = new Raycaster()
  const mouse = new Vector3()

  // Get mouse position in normalized device coordinates (-1 to +1)
  const rect = (event.target as HTMLElement).getBoundingClientRect()
  mouse.x = ((event.clientX - rect.left) / rect.width) * 2 - 1
  mouse.y = -((event.clientY - rect.top) / rect.height) * 2 + 1

  // Update raycaster
  raycaster.setFromCamera(mouse, camera)

  // Find intersections with atom meshes
  const intersections = raycaster.intersectObjects(atomMeshes, true)

  if (intersections.length > 0) {
    // Return the first intersected mesh
    const intersection = intersections[0]
    return intersection.object as THREE.Mesh
  }

  return null
}

/**
 * Convert screen coordinates to world coordinates on a plane
 * Clamps coordinates to reasonable bounds [-5, 5]
 */
export function screenToWorld(
  event: MouseEvent | React.PointerEvent,
  camera: THREE.Camera,
  planeY: number = 0
): Vector3 {
  const raycaster = new Raycaster()
  const mouse = new Vector3()
  const plane = new THREE.Plane(new THREE.Vector3(0, 1, 0), -planeY)

  // Get mouse position in normalized device coordinates
  const rect = (event.target as HTMLElement).getBoundingClientRect()
  mouse.x = ((event.clientX - rect.left) / rect.width) * 2 - 1
  mouse.y = -((event.clientY - rect.top) / rect.height) * 2 + 1

  // Update raycaster
  raycaster.setFromCamera(mouse, camera)

  // Intersect with plane
  const intersection = new Vector3()
  const hasIntersection = raycaster.ray.intersectPlane(plane, intersection)

  if (!hasIntersection) {
    // Fallback: use a fixed distance from camera
    const distance = 5
    raycaster.ray.at(distance, intersection)
  }

  // Clamp coordinates to reasonable bounds
  intersection.x = Math.max(-5, Math.min(5, intersection.x))
  intersection.y = Math.max(-5, Math.min(5, intersection.y))
  intersection.z = Math.max(-5, Math.min(5, intersection.z))

  return intersection
}

/**
 * Get world position from screen coordinates using existing ray
 */
export function getWorldPositionFromRay(
  raycaster: Raycaster,
  distance: number
): Vector3 {
  const position = new Vector3()
  raycaster.ray.at(distance, position)
  return position
}

