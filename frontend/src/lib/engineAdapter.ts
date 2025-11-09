import { MoleculeGraph } from '@biosynth/engine'
import { MoleculeSerializer } from '@biosynth/engine'
import { Atom, Bond } from '@biosynth/engine'
import { useMoleculeStore } from '../store/moleculeStore'

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

  // Update store with cloned molecule to trigger re-render
  store.setMolecule(molecule.clone())

  // Re-run predictions when geometry changes
  store.fetchPredictions()
}

/**
 * Add bond between two atoms
 */
export function addBond(a1: string, a2: string, order: number = 1): void {
  const store = useMoleculeStore.getState()
  const molecule = store.currentMolecule

  if (!molecule) return

  // Add bond
  const bondId = molecule.addBond(a1, a2, order)
  if (bondId) {
    // Update store with cloned molecule
    store.setMolecule(molecule.clone())

    // Re-run predictions when structure changes
    store.fetchPredictions()
  }
}

/**
 * Create a MoleculeGraph from renderable format
 * TODO: This function needs proper bond-to-atom mapping
 * Currently assumes bond.id contains atom IDs, which is incorrect
 */
export function renderableToMolecule(
  atoms: RenderableAtom[],
  bonds: RenderableBond[]
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
  bonds.forEach((bond) => {
    // This is incorrect - bond.id is the bond's ID, not atom IDs
    // Need to store atom IDs in the RenderableBond interface
    console.warn('renderableToMolecule: Bond mapping not implemented correctly')
  })

  return molecule
}

