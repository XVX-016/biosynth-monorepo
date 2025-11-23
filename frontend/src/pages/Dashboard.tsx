import { useEffect, useMemo, useState } from 'react'
import { motion } from 'framer-motion'
import { useNavigate } from 'react-router-dom'
import { listMolecules, convertSMILESToMolfile } from '../lib/api'
import type { MoleculeItem } from '../lib/api'
import Card from '../components/ui/Card'
import BarbellViewer from '../components/BarbellViewer'

// Caffeine SMILES: CN1C=NC2=C1C(=O)N(C(=O)N2C)C
const CAFFEINE_SMILES = 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'

/* --------- Inline Icons (clean, modern SVG) ---------- */
function SparklesIcon({ className = 'w-7 h-7 text-indigo-600' }: { className?: string }) {
  return (
    <svg className={className} viewBox="0 0 24 24" fill="none" aria-hidden="true">
      <path d="M12 2l1.5 3L17 7l-3.5 1L12 11l-1.5-3L7 7l3.5-2L12 2z" fill="currentColor" opacity="0.95"/>
      <path d="M4 14l.9 1.8L7 17l-1.1 1.2L6 20l-2-1.2L2 20l.1-1.8L1 17l2-1.2L4 14z" fill="currentColor" opacity="0.75"/>
      <path d="M20 14l.9 1.8L23 17l-1.1 1.2L22 20l-2-1.2L18 20l.1-1.8L17 17l2-1.2L20 14z" fill="currentColor" opacity="0.65"/>
    </svg>
  );
}

function CubeIcon({ className = 'w-7 h-7 text-emerald-500' }: { className?: string }) {
  return (
    <svg className={className} viewBox="0 0 24 24" fill="none" aria-hidden="true">
      <path d="M21 16V8l-9-5-9 5v8l9 5 9-5z" stroke="currentColor" strokeWidth="1.2" strokeLinecap="round" strokeLinejoin="round" fill="currentColor" opacity="0.95"/>
      <path d="M12 3v18" stroke="#ffffff" strokeWidth="0.6" strokeLinecap="round" strokeLinejoin="round" opacity="0.3"/>
      <path d="M3 8l9 5 9-5" stroke="#ffffff" strokeWidth="0.6" strokeLinecap="round" strokeLinejoin="round" opacity="0.3"/>
    </svg>
  );
}

function LibraryIcon({ className = 'w-7 h-7 text-sky-500' }: { className?: string }) {
  return (
    <svg className={className} viewBox="0 0 24 24" fill="none" aria-hidden="true">
      <path d="M3 6h18v14H3z" fill="currentColor" opacity="0.95"/>
      <path d="M7 6v14" stroke="#ffffff" strokeWidth="0.8" strokeLinecap="round" strokeLinejoin="round" opacity="0.25"/>
      <path d="M11 6v14" stroke="#ffffff" strokeWidth="0.8" strokeLinecap="round" strokeLinejoin="round" opacity="0.25"/>
      <path d="M15 6v14" stroke="#ffffff" strokeWidth="0.8" strokeLinecap="round" strokeLinejoin="round" opacity="0.25"/>
    </svg>
  );
}

export default function Dashboard(){
  const navigate = useNavigate()
  const [recent, setRecent] = useState<MoleculeItem[]>([])
  const [caffeineMolfile, setCaffeineMolfile] = useState<string | null>(null)

  // Load caffeine molfile on mount
  useEffect(() => {
    let cancelled = false
    ;(async () => {
      try {
        const result = await convertSMILESToMolfile(CAFFEINE_SMILES)
        if (!cancelled && result.molfile) {
          setCaffeineMolfile(result.molfile)
        }
      } catch (error) {
        console.error('Failed to load caffeine molfile:', error)
        if (!cancelled) setCaffeineMolfile(null)
      }
    })()
    return () => { cancelled = true }
  }, [])

  useEffect(() => {
    let cancelled = false
    ;(async () => {
      try {
        const items = await listMolecules(24)
        if (!cancelled) setRecent(items ?? [])
      } catch {
        if (!cancelled) setRecent([])
      }
    })()
    return () => { cancelled = true }
  }, [])

  const kpis = useMemo(() => {
    const total = recent.length
    const avgStability = 0.86 // placeholder
    const lastGenerated = recent[0]?.name ?? '—'
    return { total, avgStability, lastGenerated }
  }, [recent])

  return (
    <motion.div 
      initial={{ opacity: 0 }} 
      animate={{ opacity: 1 }} 
      transition={{ duration: 0.3 }} 
      className="space-y-12"
    >
      {/* Hero Section - Improved 2-column layout */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-12 items-center max-w-7xl mx-auto px-6">
        {/* Left Section - Text Content */}
        <motion.div 
          className="space-y-6"
          initial={{ opacity: 0, x: -20 }}
          animate={{ opacity: 1, x: 0 }}
          transition={{ duration: 0.6, delay: 0.1 }}
        >
          {/* Subheading */}
          <p className="text-sm font-semibold text-indigo-600 tracking-wide uppercase">
            AI-Powered Molecular Design
          </p>

          {/* Main Heading */}
          <h1 className="text-4xl md:text-5xl lg:text-6xl font-extrabold text-black leading-tight">
            Reimagine Molecular Discovery with AI
          </h1>

          {/* Description */}
          <p className="text-lg text-darkGrey leading-relaxed max-w-xl">
            Build, analyze, and explore molecular structures instantly using cutting-edge computational chemistry tools. 
            Visualize 3D structures, save to your personal library, and collaborate with the community.
          </p>

          {/* CTA Buttons */}
          <div className="flex flex-col sm:flex-row gap-4 pt-2">
            <motion.button 
              onClick={() => navigate('/lab')}
              className="btn-primary px-8 py-4 text-lg font-semibold rounded-xl shadow-lg hover:shadow-xl transition-all"
              whileHover={{ scale: 1.02 }}
              whileTap={{ scale: 0.98 }}
            >
              Generate Molecule
            </motion.button>
            <motion.button 
              onClick={() => navigate('/library')}
              className="btn-secondary px-8 py-4 text-lg font-semibold rounded-xl border-2 border-lightGrey hover:border-midGrey transition-all"
              whileHover={{ scale: 1.02 }}
              whileTap={{ scale: 0.98 }}
            >
              Explore Library
            </motion.button>
          </div>

          {/* Quick Stats */}
          <div className="flex gap-8 pt-6">
            <div>
              <div className="text-sm text-midGrey mb-1">Total saved</div>
              <div className="text-3xl font-bold text-black">{kpis.total}</div>
            </div>
            <div>
              <div className="text-sm text-midGrey mb-1">Avg stability</div>
              <div className="text-3xl font-bold text-black">{kpis.avgStability}</div>
            </div>
            <div className="hidden md:block">
              <div className="text-sm text-midGrey mb-1">Last generated</div>
              <div className="text-sm text-black truncate max-w-[120px]" title={kpis.lastGenerated}>{kpis.lastGenerated}</div>
            </div>
          </div>
        </motion.div>

        {/* Right Section - 3D Viewer (Caffeine) */}
        <motion.div
          className="flex justify-center items-center"
          initial={{ opacity: 0, scale: 0.98 }}
          animate={{ opacity: 1, scale: 1 }}
          transition={{ duration: 0.7, delay: 0.12 }}
        >
          <div className="w-full max-w-[480px] h-[420px] rounded-2xl overflow-hidden border border-lightGrey shadow-xl bg-white relative">
            {/* Subtle glow effect */}
            <div className="absolute inset-0 bg-gradient-to-br from-indigo-50/40 to-transparent pointer-events-none" />
            
            {/* 3D Viewer */}
            <div className="relative w-full h-full z-10">
              {caffeineMolfile ? (
                <BarbellViewer
                  molfile={caffeineMolfile}
                  mode="hero"
                  autorotate={true}
                  interactive={true}
                  atomScale={0.28}
                  bondRadius={0.06}
                  height={420}
                  className="w-full h-full"
                />
              ) : (
                <div className="w-full h-full flex items-center justify-center bg-gray-50">
                  <div className="text-center">
                    <div className="w-8 h-8 border-2 border-indigo-600 border-t-transparent rounded-full animate-spin mx-auto mb-2"></div>
                    <div className="text-sm text-gray-500">Loading molecule...</div>
                  </div>
                </div>
              )}
            </div>
            
            {/* Label */}
            <div className="absolute bottom-4 left-1/2 transform -translate-x-1/2 bg-black/75 text-white px-4 py-2 rounded-full text-sm backdrop-blur-sm">
              Molecule Preview: Caffeine (C₈H₁₀N₄O₂)
            </div>
          </div>
        </motion.div>
      </div>
      
      {/* Features Section */}
      <section className="max-w-7xl mx-auto px-6">
        <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
          <Card className="p-6 hover:shadow-neon-hover transition-all duration-300 flex gap-4 items-start">
            <div className="mt-1 flex-shrink-0">
              <SparklesIcon />
            </div>
            <div>
              <h3 className="text-xl font-semibold text-black mb-2">AI Molecule Generator</h3>
              <p className="text-darkGrey text-sm leading-relaxed">
                Build molecules instantly using machine learning. Generate novel structures from natural language prompts.
              </p>
            </div>
          </Card>
          <Card className="p-6 hover:shadow-neon-hover transition-all duration-300 flex gap-4 items-start">
            <div className="mt-1 flex-shrink-0">
              <CubeIcon />
            </div>
            <div>
              <h3 className="text-xl font-semibold text-black mb-2">3D Molecular Viewer</h3>
              <p className="text-darkGrey text-sm leading-relaxed">
                Interact with high-quality 3D visualizations. Rotate, zoom, and explore molecular structures in real-time.
              </p>
            </div>
          </Card>
          <Card className="p-6 hover:shadow-neon-hover transition-all duration-300 flex gap-4 items-start">
            <div className="mt-1 flex-shrink-0">
              <LibraryIcon />
            </div>
            <div>
              <h3 className="text-xl font-semibold text-black mb-2">Your Molecule Library</h3>
              <p className="text-darkGrey text-sm leading-relaxed">
                Save, manage, and explore your creations. Build a personal collection of molecular structures.
              </p>
            </div>
          </Card>
        </div>
      </section>

      {/* Stats Grid - Compact */}
      <section className="max-w-7xl mx-auto px-6">
        <div className="grid grid-cols-1 sm:grid-cols-3 gap-6">
          <Card className="p-6">
            <div className="text-sm text-midGrey mb-2">Total saved</div>
            <div className="text-4xl font-bold text-black">{kpis.total}</div>
            <div className="text-xs text-midGrey mt-2">Molecules in library</div>
          </Card>
          <Card className="p-6">
            <div className="text-sm text-midGrey mb-2">Avg stability</div>
            <div className="text-4xl font-bold text-black">{kpis.avgStability}</div>
            <div className="mt-3 h-6">
              <svg viewBox="0 0 100 24" className="w-full h-full text-darkGrey">
                <polyline fill="none" stroke="currentColor" strokeWidth="2" points="0,20 10,18 20,12 30,14 40,9 50,10 60,7 70,9 80,6 90,8 100,4" />
              </svg>
            </div>
          </Card>
          <Card className="p-6">
            <div className="text-sm text-midGrey mb-2">Last generated</div>
            <div className="text-2xl font-bold truncate text-black" title={kpis.lastGenerated}>{kpis.lastGenerated}</div>
            <div className="text-xs text-midGrey mt-2">Most recent creation</div>
          </Card>
        </div>
      </section>

      {/* Recent Molecules - Horizontal Scroll */}
      <section className="max-w-7xl mx-auto px-6">
        <div className="flex items-center justify-between mb-6">
          <h2 className="text-2xl font-semibold text-black">Recent Molecules</h2>
          <button 
            onClick={() => navigate('/library')}
            className="text-sm text-midGrey hover:text-black transition-colors"
          >
            View all →
          </button>
        </div>
        
        <div className="flex overflow-x-auto gap-6 pb-4 -mx-6 px-6 scrollbar-hide">
          {recent.length > 0 ? (
            recent.slice(0, 12).map((m) => (
              <motion.div
                key={m.id}
                className="min-w-[220px] flex-shrink-0"
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ duration: 0.28 }}
              >
                <div 
                  className="cursor-pointer h-full"
                  onClick={() => navigate(`/lab?id=${m.id}`)}
                >
                <Card className="overflow-hidden hover:shadow-neon-hover transition-all h-full">
                  <div className="h-36 flex items-center justify-center overflow-hidden bg-offwhite">
                    {m.molfile ? (
                      <BarbellViewer
                        molfile={m.molfile}
                        mode="card"
                        height={140}
                        atomScale={0.18}
                        bondRadius={0.04}
                      />
                    ) : m.thumbnail_b64 ? (
                      <img 
                        src={m.thumbnail_b64.startsWith('data:') ? m.thumbnail_b64 : `data:image/png;base64,${m.thumbnail_b64}`} 
                        alt={m.name} 
                        className="max-h-full max-w-full object-contain" 
                      />
                    ) : (
                      <div className="text-midGrey text-sm">No preview</div>
                    )}
                  </div>
                  <div className="px-4 py-3">
                    <div className="font-medium truncate text-black">{m.name}</div>
                    {m.formula && (
                      <div className="text-xs text-midGrey mt-1 font-mono">{m.formula}</div>
                    )}
                    <div className="text-xs text-midGrey mt-1">{new Date(m.created_at).toLocaleDateString()}</div>
                  </div>
                </Card>
                </div>
              </motion.div>
            ))
          ) : (
            <div className="text-midGrey col-span-full text-center py-12 w-full">
              <p className="mb-4">No molecules yet.</p>
              <button 
                onClick={() => navigate('/lab')}
                className="btn-primary px-6 py-2"
              >
                Generate Your First Molecule
              </button>
            </div>
          )}
        </div>
      </section>

      {/* AI Model Cards */}
      <section className="max-w-7xl mx-auto px-6">
        <h2 className="text-2xl font-semibold text-black mb-6">AI Models</h2>
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-6">
          <Card className="p-6 hover:shadow-neon-hover transition-shadow">
            <h3 className="text-lg font-semibold text-black mb-2">CoreGen-1</h3>
            <p className="text-sm text-darkGrey">Molecule generation from SMILES</p>
          </Card>
          <Card className="p-6 hover:shadow-neon-hover transition-shadow">
            <h3 className="text-lg font-semibold text-black mb-2">PropNet-X</h3>
            <p className="text-sm text-darkGrey">Property prediction</p>
          </Card>
          <Card className="p-6 hover:shadow-neon-hover transition-shadow">
            <h3 className="text-lg font-semibold text-black mb-2">ReactFlow-R1</h3>
            <p className="text-sm text-darkGrey">Synthesis planning</p>
          </Card>
          <Card className="p-6 hover:shadow-neon-hover transition-shadow">
            <h3 className="text-lg font-semibold text-black mb-2">ChemGPT-S</h3>
            <p className="text-sm text-darkGrey">Reaction simulation</p>
          </Card>
        </div>
      </section>
    </motion.div>
  )
}
