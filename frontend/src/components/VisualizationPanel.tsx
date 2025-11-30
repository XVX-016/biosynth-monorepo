import { useEffect } from 'react'
import { useMoleculeStore, RenderMode, ColorScheme } from '../store/moleculeStore'

const renderModes: { id: RenderMode; label: string }[] = [
  { id: 'ballstick', label: 'Ball & Stick' },
  { id: 'spacefill', label: 'CPK' },
  { id: 'wireframe', label: 'Wireframe' },
]

const colorSchemes: { id: ColorScheme; label: string }[] = [
  { id: 'element', label: 'Element' },
  { id: 'monochrome', label: 'Monochrome' },
  { id: 'electrostatic', label: 'Electrostatic' },
]

export default function VisualizationPanel() {
  const renderMode = useMoleculeStore((s) => s.renderMode)
  const colorScheme = useMoleculeStore((s) => s.colorScheme)
  const setRenderMode = useMoleculeStore((s) => s.setRenderMode)
  const setColorScheme = useMoleculeStore((s) => s.setColorScheme)

  useEffect(() => {
    if (typeof window === 'undefined') return
    const storedMode = window.localStorage.getItem('molforge-render-mode') as RenderMode | null
    const storedScheme = window.localStorage.getItem('molforge-color-scheme') as ColorScheme | null
    if (storedMode && storedMode !== renderMode) {
      setRenderMode(storedMode)
    }
    if (storedScheme && storedScheme !== colorScheme) {
      setColorScheme(storedScheme)
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])

  useEffect(() => {
    if (typeof window === 'undefined') return
    window.localStorage.setItem('molforge-render-mode', renderMode)
  }, [renderMode])

  useEffect(() => {
    if (typeof window === 'undefined') return
    window.localStorage.setItem('molforge-color-scheme', colorScheme)
  }, [colorScheme])

  return (
    <div className="space-y-4">
      <div className="flex items-center justify-between">
        <h3 className="text-base font-semibold text-black">Visualization</h3>
        <span className="text-xs uppercase tracking-[0.3em] text-midGrey">Modes</span>
      </div>
      <div className="space-y-2">
        <p className="text-xs uppercase tracking-[0.2em] text-midGrey">Model</p>
        <div className="grid grid-cols-3 gap-2">
          {renderModes.map((mode) => (
            <button
              key={mode.id}
              className={`text-xs font-semibold rounded-lg border px-2 py-2 transition ${
                renderMode === mode.id
                  ? 'bg-black text-white border-black shadow-neon-sm'
                  : 'bg-white text-darkGrey border-lightGrey hover:border-darkGrey/40'
              }`}
              onClick={() => setRenderMode(mode.id)}
            >
              {mode.label}
            </button>
          ))}
        </div>
      </div>
      <div className="space-y-2">
        <p className="text-xs uppercase tracking-[0.2em] text-midGrey">Color</p>
        <div className="grid grid-cols-3 gap-2">
          {colorSchemes.map((scheme) => (
            <button
              key={scheme.id}
              className={`text-xs font-semibold rounded-lg border px-2 py-2 transition ${
                colorScheme === scheme.id
                  ? 'bg-black text-white border-black shadow-neon-sm'
                  : 'bg-white text-darkGrey border-lightGrey hover:border-darkGrey/40'
              }`}
              onClick={() => setColorScheme(scheme.id)}
            >
              {scheme.label}
            </button>
          ))}
        </div>
      </div>
    </div>
  )
}



