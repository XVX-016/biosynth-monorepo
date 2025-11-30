import type { MoleculeGraph } from '@biosynth/engine'
import { MoleculeSerializer } from '@biosynth/engine'
import { apiClient } from '../lib/api'

interface MlAtomPayload {
  id: string
  element: string
  position: [number, number, number]
}

interface MlBondPayload {
  id: string
  atom1: string
  atom2: string
  order: number
}

export interface GnassResponse {
  predictions?: Record<string, number>
  explanation?: string
  features?: Record<string, number>
}

const mapMoleculeToPayload = (molecule: MoleculeGraph) => {
  const atoms: MlAtomPayload[] = []
  const bonds: MlBondPayload[] = []

  molecule.atoms.forEach((atom) => {
    atoms.push({
      id: atom.id,
      element: atom.element,
      position: atom.position,
    })
  })

  molecule.bonds.forEach((bond) => {
    bonds.push({
      id: bond.id,
      atom1: bond.a1,
      atom2: bond.a2,
      order: bond.order,
    })
  })

  return { atoms, bonds }
}

export async function runGnassPipeline(molecule: MoleculeGraph): Promise<GnassResponse> {
  const payload = mapMoleculeToPayload(molecule)
  const response = await apiClient.post<GnassResponse>('/api/ml/predict', payload)
  return response.data
}

export async function runGnassWithEngines(
  molecule: MoleculeGraph,
  engineOutputs?: Record<string, unknown>
): Promise<GnassResponse> {
  const payload = {
    molecule: mapMoleculeToPayload(molecule),
    engine_outputs: engineOutputs,
    use_engine_features: Boolean(engineOutputs),
  }

  const response = await apiClient.post<GnassResponse>('/api/ml/predict/with-engines', payload)
  return response.data
}

export function moleculeFingerprint(molecule: MoleculeGraph | null): string {
  if (!molecule) return 'empty'
  try {
    const json = MoleculeSerializer.toJSON(molecule)
    return JSON.stringify({ atoms: json.atoms.length, bonds: json.bonds.length, checksum: json.atoms.slice(0, 6) })
  } catch {
    return `${molecule.atoms.size}-${molecule.bonds.size}`
  }
}


