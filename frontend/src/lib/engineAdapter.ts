import { MoleculeGraph, MoleculeSerializer, ForceField, autoBondNewAtom } from '@biosynth/engine'
import { useMoleculeStore } from '../store/moleculeStore'
import { pushState } from '../store/historyStore'

export interface RenderableAtom {
  id: string
  element: string
  position: [number, number, number]
}

export interface RenderableBond {
  id: string
  from: [number, number, number]
  to: [number, number, number]
  order: number
}

/**
 * Converts MoleculeGraph to React-renderable format
 */
export function moleculeToRenderable(molecule: MoleculeGraph | null): {
  atoms: RenderableAtom[]
  bonds: RenderableBond[]
  smiles: string
} {
  if (!molecule) {
    return { atoms: [], bonds: [], smiles: '' }
  }

  const atoms: RenderableAtom[] = Array.from(molecule.atoms.values()).map(
    (atom) => ({
      id: atom.id,
      element: atom.element,
      position: atom.position,
    })
  )

  const bonds: RenderableBond[] = Array.from(molecule.bonds.values()).map(
    (bond) => {
      const atom1 = molecule.atoms.get(bond.a1)
      const atom2 = molecule.atoms.get(bond.a2)

      if (!atom1 || !atom2) {
        throw new Error(`Bond references missing atoms: ${bond.a1}, ${bond.a2}`)
      }

      return {
        id: bond.id,
        from: atom1.position,
        to: atom2.position,
        order: bond.order,
      }
    }
  )

  // Convert to SMILES when structure changes
  const smiles = MoleculeSerializer.toSMILES(molecule)

  return { atoms, bonds, smiles }
}

/**
 * Update atom position in molecule graph
 */
export function updateAtomPosition(
  id: string,
  position: [number, number, number]
): void {
  const store = useMoleculeStore.getState()
  const molecule = store.currentMolecule

  if (!molecule) return

  const atom = molecule.atoms.get(id)
  if (!atom) return

  // Update position
  molecule.atoms.set(id, { ...atom, position })

  // TODO: Run geometry optimization in WebWorker for better performance
  // For now, run synchronously with a small number of iterations
  // This will relax the molecule after dragging
  setTimeout(() => {
    const cloned = molecule.clone()
    ForceField.optimizeGeometry(cloned, 5, 0.005) // Quick relaxation
    store.setMolecule(cloned)
  }, 0)

  // Update store with cloned molecule to trigger re-render
  store.setMolecule(molecule.clone())

  // Re-run predictions when geometry changes
  store.fetchPredictions()
}

/**
 * Add atom to molecule graph
 */
export function addAtom(
  element: string,
  position: [number, number, number]
): void {
  const store = useMoleculeStore.getState()
  const molecule = store.currentMolecule

  if (!molecule) {
    // Create new molecule if none exists
    const newMolecule = new MoleculeGraph()
    newMolecule.addAtom({ element: element as any, position })
    store.setMolecule(newMolecule)
    pushState()
    return
  }

  // Clone first to avoid mutation
  const cloned = molecule.clone()

  // Add atom to cloned molecule
  const newId = cloned.addAtom({ element: element as any, position })

  // Auto-bond if enabled and we have more than 1 atom
  if (store.autoBond && newId && cloned.atoms.size > 1) {
    console.log('[Auto-Bond] Attempting to bond new atom:', newId)
    autoBondNewAtom(cloned, newId)
    console.log('[Auto-Bond] Bonds after auto-bond:', cloned.bonds.size)
  }

  // Update store with bonded molecule
  store.setMolecule(cloned)
  pushState()

  // Optimize geometry after a short delay
  setTimeout(() => {
    const optimized = cloned.clone()
    ForceField.optimizeGeometry(optimized, 10, 0.01)
    store.setMolecule(optimized)
    // Re-run predictions after optimization
    store.fetchPredictions()
  }, 50)
}

/**
 * Remove atom from molecule graph
 */
export function removeAtom(id: string): void {
  const store = useMoleculeStore.getState()
  const molecule = store.currentMolecule

  if (!molecule) return

  pushState()
  molecule.removeAtom(id)
  const cloned = molecule.clone()
  store.setMolecule(cloned)

  // Clear selection if deleted atom was selected
  if (store.selectedAtomId === id) {
    store.selectAtom(null)
  }

  // Re-run predictions
  store.fetchPredictions()
}

/**
 * Remove bond from molecule graph
 */
export function removeBond(id: string): void {
  const store = useMoleculeStore.getState()
  const molecule = store.currentMolecule

  if (!molecule) return

  pushState()
  molecule.removeBond(id)
  const cloned = molecule.clone()
  store.setMolecule(cloned)

  // Clear selection if deleted bond was selected
  if (store.selectedBondId === id) {
    store.selectBond(null)
  }

  // Re-run predictions
  store.fetchPredictions()
}

/**
 * Add bond between two atoms
 */
export function addBond(a1: string, a2: string, order: number = 1): void {
  const store = useMoleculeStore.getState()
  const molecule = store.currentMolecule

  if (!molecule) return

  pushState()

  // Add bond
  const bondId = molecule.addBond(a1, a2, order)
  if (bondId) {
    // Update store with cloned molecule
    const cloned = molecule.clone()
    store.setMolecule(cloned)

    // TODO: Run geometry optimization in WebWorker for better performance
    // Optimize geometry after bond addition to relax the structure
    setTimeout(() => {
      const optimized = cloned.clone()
      ForceField.optimizeGeometry(optimized, 10, 0.005) // More iterations for new bonds
      store.setMolecule(optimized)
    }, 0)

    // Re-run predictions when structure changes
    store.fetchPredictions()
  }
}


/**
 * Convert MoleculeGraph to JSON for saving
 */
export function moleculeToJSON(molecule: MoleculeGraph | null): string {
  if (!molecule) return JSON.stringify({ atoms: [], bonds: [] })
  return JSON.stringify(MoleculeSerializer.toJSON(molecule))
}

/**
 * Create MoleculeGraph from JSON
 */
export function moleculeFromJSON(jsonStr: string): MoleculeGraph | null {
  try {
    const data = JSON.parse(jsonStr)
    return MoleculeSerializer.fromJSON(data)
  } catch (error) {
    console.error('Failed to parse molecule JSON:', error)
    return null
  }
}

/**
 * Get thumbnail from canvas (base64)
 */
export function getCanvasThumbnail(): string | null {
  const canvas = document.querySelector('canvas')
  if (!canvas) return null
  return canvas.toDataURL('image/png')
}

/**
 * Create a MoleculeGraph from renderable format
 * TODO: This function needs proper bond-to-atom mapping
 * Currently assumes bond.id contains atom IDs, which is incorrect
 */
export function renderableToMolecule(
  atoms: RenderableAtom[],
  _bonds: RenderableBond[] // Unused until bond mapping is fixed
): MoleculeGraph {
  const molecule = new MoleculeGraph()
  const atomIdMap = new Map<string, string>() // old ID -> new ID

  // Add atoms and track ID mapping
  atoms.forEach((atom) => {
    const newId = molecule.addAtom({
      element: atom.element as any,
      position: atom.position,
    })
    atomIdMap.set(atom.id, newId)
  })

  // Add bonds using mapped IDs
  // TODO: Fix bond mapping - bonds need to reference atoms by position or store atom IDs
  // For now, this is a placeholder that won't work correctly
  // Note: bonds parameter is intentionally unused until bond mapping is fixed

  return molecule
}

