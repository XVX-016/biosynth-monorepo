export type SpectroscopyStatus = 'idle' | 'calculating' | 'complete' | 'error'

export interface IRSpectrumPeak {
  frequency: number
  intensity: number
  bondType: string
  atoms: string[]
  annotation?: string
}

export interface IRSpectrum {
  peaks: IRSpectrumPeak[]
  spectrumImage?: string
}

export interface ProtonNMRPeak {
  chemicalShift: number
  integration: number
  multiplicity: string
  atoms: string[]
  annotation?: string
}

export interface CarbonNMRPeak {
  chemicalShift: number
  atoms: string[]
  annotation?: string
}

export interface NMRSpectrum {
  protons: ProtonNMRPeak[]
  carbon13?: CarbonNMRPeak[]
}

export interface MassSpectrumPeak {
  mz: number
  intensity: number
  formula: string
  fragment: string
  atoms?: string[]
}

export interface MassSpectrum {
  peaks: MassSpectrumPeak[]
  parentMass: number
}

export interface SpectralData {
  ir: IRSpectrum | null
  nmr: NMRSpectrum | null
  mass: MassSpectrum | null
  status: SpectroscopyStatus
  error?: string | null
  timestamp?: string
}



