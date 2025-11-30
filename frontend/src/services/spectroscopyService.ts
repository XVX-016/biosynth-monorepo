import type { MoleculeGraph } from '@biosynth/engine'
import type {
  IRSpectrum,
  NMRSpectrum,
  MassSpectrum,
  IRSpectrumPeak,
  ProtonNMRPeak,
  CarbonNMRPeak,
  MassSpectrumPeak,
} from '../types/spectroscopy'
import {
  IR_FREQUENCY_WINDOWS,
  NMR_SHIFT_REFERENCES,
  getAtomicMass,
} from '../utils/chemistry'

interface SimpleAtom {
  id: string
  element: string
  position: [number, number, number]
  neighbors: string[]
}

interface SimpleBond {
  id: string
  a1: string
  a2: string
  order: number
}

export class SpectroscopyService {
  async calculateAll(molecule: MoleculeGraph): Promise<{
    ir: IRSpectrum
    nmr: NMRSpectrum
    mass: MassSpectrum
  }> {
    const simple = this.toSimpleMolecule(molecule)
    return {
      ir: this.calculateIR(simple),
      nmr: this.calculateNMR(simple),
      mass: this.calculateMass(simple),
    }
  }

  private toSimpleMolecule(molecule: MoleculeGraph) {
    const atoms: Record<string, SimpleAtom> = {}
    molecule.atoms.forEach((atom, id) => {
      atoms[id] = {
        id,
        element: atom.element,
        position: atom.position,
        neighbors: [],
      }
    })
    const bonds: SimpleBond[] = []
    molecule.bonds.forEach((bond, id) => {
      bonds.push({ id, a1: bond.a1, a2: bond.a2, order: bond.order })
      atoms[bond.a1]?.neighbors.push(bond.a2)
      atoms[bond.a2]?.neighbors.push(bond.a1)
    })
    return { atoms: Object.values(atoms), bonds }
  }

  private calculateIR(simple: { atoms: SimpleAtom[]; bonds: SimpleBond[] }): IRSpectrum {
    const peaks: IRSpectrumPeak[] = []

    simple.bonds.forEach((bond) => {
      const atom1 = simple.atoms.find((a) => a.id === bond.a1)
      const atom2 = simple.atoms.find((a) => a.id === bond.a2)
      if (!atom1 || !atom2) return

      const pair = [atom1.element, atom2.element].sort().join('-')
      const orderSymbol = bond.order === 3 ? '≡' : bond.order === 2 ? '=' : '-'
      const descriptor =
        pair === 'C-O' && bond.order === 2
          ? 'C=O'
          : pair === 'C-O' && bond.order === 1
          ? 'C-O'
          : pair === 'C-N' && bond.order === 2
          ? 'C=N'
          : `${pair.replace('-', orderSymbol)}`

      const window = IR_FREQUENCY_WINDOWS[descriptor] ?? IR_FREQUENCY_WINDOWS[pair]
      if (!window) return

      const frequency =
        window.min + Math.random() * (window.max - window.min)
      const intensity = 0.6 + Math.random() * 0.4

      peaks.push({
        frequency: parseFloat(frequency.toFixed(0)),
        intensity: parseFloat(intensity.toFixed(2)),
        bondType: descriptor,
        atoms: [atom1.id, atom2.id],
        annotation: window.label,
      })
    })

    // Add broad O-H / N-H if present
    const hasOH = simple.bonds.some((bond) => {
      const a1 = simple.atoms.find((a) => a.id === bond.a1)
      const a2 = simple.atoms.find((a) => a.id === bond.a2)
      return (
        bond.order === 1 &&
        ((a1?.element === 'O' && a2?.element === 'H') ||
          (a1?.element === 'N' && a2?.element === 'H'))
      )
    })
    if (hasOH) {
      peaks.push({
        frequency: 3400 + Math.random() * 200,
        intensity: 0.8,
        bondType: 'O-H',
        atoms: [],
        annotation: 'Hydrogen-bonded stretch',
      })
    }

    return { peaks: peaks.sort((a, b) => b.frequency - a.frequency) }
  }

  private calculateNMR(simple: { atoms: SimpleAtom[]; bonds: SimpleBond[] }): NMRSpectrum {
    const protons: ProtonNMRPeak[] = []
    const carbons: CarbonNMRPeak[] = []

    simple.atoms.forEach((atom) => {
      if (atom.element === 'H') {
        const parent = simple.atoms.find((a) => a.id === atom.neighbors[0])
        if (!parent) return
        const env = this.inferProtonEnvironment(parent)
        const range = NMR_SHIFT_REFERENCES[env] ?? { min: 1, max: 5 }
        const shift =
          range.min + Math.random() * (range.max - range.min)
        protons.push({
          chemicalShift: parseFloat(shift.toFixed(2)),
          integration: 1,
          multiplicity: this.estimateMultiplicity(parent.neighbors.length),
          atoms: [atom.id],
          annotation: env,
        })
      }

      if (atom.element === 'C') {
        const neighborElements = atom.neighbors
          .map((id) => simple.atoms.find((a) => a.id === id)?.element)
          .filter(Boolean) as string[]
        const env = neighborElements.includes('O')
          ? 'carbonyl'
          : neighborElements.includes('N')
          ? 'amine'
          : neighborElements.filter((el) => el === 'C').length >= 2
          ? 'aromatic'
          : 'alkyl'
        const base =
          env === 'carbonyl'
            ? 170 + Math.random() * 20
            : env === 'aromatic'
            ? 120 + Math.random() * 20
            : env === 'amine'
            ? 40 + Math.random() * 20
            : 20 + Math.random() * 40
        carbons.push({
          chemicalShift: parseFloat(base.toFixed(1)),
          atoms: [atom.id],
          annotation: env,
        })
      }
    })

    return {
      protons: protons.sort((a, b) => b.chemicalShift - a.chemicalShift),
      carbon13: carbons.sort((a, b) => b.chemicalShift - a.chemicalShift),
    }
  }

  private calculateMass(simple: { atoms: SimpleAtom[]; bonds: SimpleBond[] }): MassSpectrum {
    const parentMass = simple.atoms.reduce(
      (acc, atom) => acc + getAtomicMass(atom.element),
      0
    )

    const peaks: MassSpectrumPeak[] = [
      {
        mz: parseFloat(parentMass.toFixed(1)),
        intensity: 1,
        formula: 'M⁺',
        fragment: 'parent',
      },
    ]

    // Generate a few fragment peaks
    const fragmentCount = Math.min(4, Math.max(1, simple.atoms.length / 5))
    for (let i = 0; i < fragmentCount; i++) {
      const retainedAtoms = simple.atoms.filter(
        () => Math.random() > 0.3
      )
      if (retainedAtoms.length === 0) continue
      const mass = retainedAtoms.reduce(
        (acc, atom) => acc + getAtomicMass(atom.element),
        0
      )
      peaks.push({
        mz: parseFloat(mass.toFixed(1)),
        intensity: 0.3 + Math.random() * 0.6,
        formula: `M-${i + 1}`,
        fragment: 'fragment',
        atoms: retainedAtoms.map((atom) => atom.id),
      })
    }

    return {
      peaks: peaks
        .filter((peak, index, arr) => {
          // remove duplicates by mz
          return (
            arr.findIndex((p) => Math.abs(p.mz - peak.mz) < 0.1) === index
          )
        })
        .sort((a, b) => a.mz - b.mz),
      parentMass: parseFloat(parentMass.toFixed(1)),
    }
  }

  private inferProtonEnvironment(atom: SimpleAtom): keyof typeof NMR_SHIFT_REFERENCES {
    if (atom.element === 'O' || atom.element === 'N') return 'alcohol'
    const carbonNeighbors = atom.neighbors.length
    if (carbonNeighbors >= 3) return 'aromatic'
    if (carbonNeighbors === 2) return 'allylic'
    return 'alkyl'
  }

  private estimateMultiplicity(neighborCount: number): string {
    const translation: Record<number, string> = {
      0: 's',
      1: 'd',
      2: 't',
      3: 'q',
    }
    return translation[Math.min(3, neighborCount)] ?? 'm'
  }
}


