import { useState } from 'react'
import { MoleculeSerializer } from '@biosynth/engine'
import { useMoleculeStore } from '../store/moleculeStore'
import { getCanvasThumbnail, moleculeToJSON } from '../lib/engineAdapter'
import { downloadDataUrl, downloadText } from '../utils/download'
import { convertSMILESToMolfile } from '../lib/api'
import { encodeSharePayload } from '../utils/shareLink'

type ExportState = 'idle' | 'mol'

export default function ExportPanel() {
  const molecule = useMoleculeStore((s) => s.currentMolecule)
  const [status, setStatus] = useState<ExportState>('idle')

  const ensureMolecule = () => {
    if (!molecule) {
      alert('Add atoms to enable exporting.')
      return false
    }
    return true
  }

  const handleCopySmiles = async () => {
    if (!ensureMolecule()) return
    const smiles = MoleculeSerializer.toSMILES(molecule)
    if (!smiles) {
      alert('Unable to serialize molecule to SMILES.')
      return
    }
    await navigator.clipboard.writeText(smiles)
    alert('SMILES copied to clipboard.')
  }

  const handleExportPNG = () => {
    if (!ensureMolecule()) return
    const dataUrl = getCanvasThumbnail()
    if (!dataUrl) {
      alert('Canvas not ready. Move the molecule once and try again.')
      return
    }
    downloadDataUrl('molforge-lab.png', dataUrl)
  }

  const handleExportSVG = () => {
    if (!ensureMolecule()) return
    const dataUrl = getCanvasThumbnail()
    if (!dataUrl) {
      alert('Canvas not ready. Move the molecule once and try again.')
      return
    }
    const svg = `<svg xmlns="http://www.w3.org/2000/svg" width="1280" height="720" viewBox="0 0 1280 720">
  <rect width="1280" height="720" fill="#ffffff"/>
  <image href="${dataUrl}" x="0" y="0" width="1280" height="720" preserveAspectRatio="xMidYMid meet"/>
</svg>`
    downloadText('molforge-lab.svg', svg, 'image/svg+xml')
  }

  const handleExportMolfile = async () => {
    if (!ensureMolecule()) return
    const smiles = MoleculeSerializer.toSMILES(molecule)
    if (!smiles) {
      alert('Unable to serialize molecule to SMILES.')
      return
    }

    setStatus('mol')
    try {
      const { molfile } = await convertSMILESToMolfile(smiles)
      downloadText('molforge-export.mol', molfile, 'chemical/x-mdl-molfile')
    } catch (error) {
      console.error(error)
      alert('Failed to generate molfile. Ensure the backend is running.')
    } finally {
      setStatus('idle')
    }
  }

  const handleShareLink = async () => {
    if (!ensureMolecule()) return
    if (typeof window === 'undefined') return
    const json = moleculeToJSON(molecule)
    const encoded = encodeSharePayload(json)
    const url = `${window.location.origin}/lab?share=${encoded}`
    await navigator.clipboard.writeText(url)
    alert('Shareable link copied to clipboard.')
  }

  return (
    <div className="space-y-4">
      <div className="flex items-center justify-between">
        <h3 className="text-base font-semibold text-black">Export & Share</h3>
        <p className="text-xs uppercase tracking-[0.3em] text-midGrey">Phase 1</p>
      </div>
      <div className="grid grid-cols-2 gap-2">
        <button className="btn-secondary text-sm py-2" onClick={handleCopySmiles}>
          Copy SMILES
        </button>
        <button className="btn-secondary text-sm py-2" onClick={handleShareLink}>
          Share Link
        </button>
      </div>
      <div className="grid grid-cols-2 gap-2">
        <button className="btn-secondary text-sm py-2" onClick={handleExportPNG}>
          Export PNG
        </button>
        <button className="btn-secondary text-sm py-2" onClick={handleExportSVG}>
          Export SVG
        </button>
      </div>
      <button
        className="btn-secondary w-full text-sm py-2 disabled:opacity-50"
        onClick={handleExportMolfile}
        disabled={status === 'mol'}
      >
        {status === 'mol' ? 'Generating Molfile…' : 'Download MOL'}
      </button>
    </div>
  )
}



