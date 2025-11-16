// src/templates/placement.ts
import { MoleculeGraph } from "@biosynth/engine";
import { getBondLength } from "@biosynth/engine";

/**
 * Suggest a bond vector (direction) for attaching a new group to an atom.
 * This estimates the best direction to place new atoms to avoid overlaps.
 */
export function suggestBondVector(
  graph: MoleculeGraph,
  attachAtomId: string
): { x: number; y: number; z: number } {
  const attachAtom = graph.atoms.get(attachAtomId);
  if (!attachAtom) {
    // Default direction if atom not found
    return { x: 1, y: 0, z: 0 };
  }

  // Get all neighbors of the attachment atom
  const neighbors = graph.getNeighbors(attachAtomId);
  
  if (neighbors.length === 0) {
    // No neighbors - use default direction
    return { x: 1, y: 0, z: 0 };
  }

  // Calculate average direction away from neighbors
  let sumX = 0;
  let sumY = 0;
  let sumZ = 0;

  for (const neighborId of neighbors) {
    const neighbor = graph.atoms.get(neighborId);
    if (!neighbor) continue;

    // Vector from attachment atom to neighbor
    const dx = neighbor.position[0] - attachAtom.position[0];
    const dy = neighbor.position[1] - attachAtom.position[1];
    const dz = neighbor.position[2] - attachAtom.position[2];

    // Normalize and accumulate (we want opposite direction)
    const len = Math.sqrt(dx * dx + dy * dy + dz * dz);
    if (len > 0.01) {
      sumX -= dx / len;
      sumY -= dy / len;
      sumZ -= dz / len;
    }
  }

  // Normalize the result
  const len = Math.sqrt(sumX * sumX + sumY * sumY + sumZ * sumZ);
  if (len > 0.01) {
    return {
      x: sumX / len,
      y: sumY / len,
      z: sumZ / len,
    };
  }

  // Fallback: use perpendicular direction
  // If neighbors are in a plane, try to find a perpendicular vector
  if (neighbors.length >= 2) {
    const n1 = graph.atoms.get(neighbors[0])!;
    const n2 = graph.atoms.get(neighbors[1])!;
    
    const v1 = [
      n1.position[0] - attachAtom.position[0],
      n1.position[1] - attachAtom.position[1],
      n1.position[2] - attachAtom.position[2],
    ];
    const v2 = [
      n2.position[0] - attachAtom.position[0],
      n2.position[1] - attachAtom.position[1],
      n2.position[2] - attachAtom.position[2],
    ];

    // Cross product to get perpendicular direction
    const crossX = v1[1] * v2[2] - v1[2] * v2[1];
    const crossY = v1[2] * v2[0] - v1[0] * v2[2];
    const crossZ = v1[0] * v2[1] - v1[1] * v2[0];

    const crossLen = Math.sqrt(crossX * crossX + crossY * crossY + crossZ * crossZ);
    if (crossLen > 0.01) {
      return {
        x: crossX / crossLen,
        y: crossY / crossLen,
        z: crossZ / crossLen,
      };
    }
  }

  // Final fallback
  return { x: 1, y: 0, z: 0 };
}

/**
 * Place template atoms relative to an attachment point.
 * Uses suggestBondVector to determine optimal placement direction.
 */
export function placeTemplate(
  graph: MoleculeGraph,
  newAtomIds: string[],
  attachAtomId: string
): void {
  const attachAtom = graph.atoms.get(attachAtomId);
  if (!attachAtom) return;

  const direction = suggestBondVector(graph, attachAtomId);
  const bondLength = getBondLength(attachAtom.element, "C"); // Default to C-C bond length

  // Place first atom along the suggested direction
  if (newAtomIds.length > 0) {
    const firstAtom = graph.atoms.get(newAtomIds[0]);
    if (firstAtom) {
      const distance = bondLength + 0.3; // Add small spacing
      firstAtom.position[0] = attachAtom.position[0] + direction.x * distance;
      firstAtom.position[1] = attachAtom.position[1] + direction.y * distance;
      firstAtom.position[2] = attachAtom.position[2] + direction.z * distance;
    }

    // Place subsequent atoms relative to the first
    for (let i = 1; i < newAtomIds.length; i++) {
      const prevAtom = graph.atoms.get(newAtomIds[i - 1]);
      const currAtom = graph.atoms.get(newAtomIds[i]);
      if (!prevAtom || !currAtom) continue;

      // Use a simple offset (can be improved with proper geometry)
      const offset = 1.2; // Standard bond distance
      currAtom.position[0] = prevAtom.position[0] + direction.x * offset * 0.5;
      currAtom.position[1] = prevAtom.position[1] + direction.y * offset * 0.5;
      currAtom.position[2] = prevAtom.position[2] + direction.z * offset * 0.5;
    }
  }
}

