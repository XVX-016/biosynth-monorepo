import { useEffect, useMemo, useState } from 'react'
import { motion } from 'framer-motion'
import { useNavigate } from 'react-router-dom'
import BenzeneGLBViewer from '../components/BenzeneGLBViewer'
import { listMolecules } from '../lib/api'
import type { MoleculeItem } from '../lib/api'
import Card from '../components/ui/Card'
import BarbellViewer from '../components/BarbellViewer'

export default function Dashboard(){
  const navigate = useNavigate()
  const [recent, setRecent] = useState<MoleculeItem[]>([])

  useEffect(() => {
    let cancelled = false
    ;(async () => {
      try {
        const items = await listMolecules(12)
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
          <p className="text-sm font-semibold text-blue-600 tracking-wide uppercase">
            AI-Powered Molecular Design
          </p>

          {/* Main Heading */}
          <h1 className="text-4xl md:text-5xl lg:text-6xl font-bold text-black leading-tight">
            Reimagine Molecular<br />Discovery with AI
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
          <div className="flex gap-8 pt-4">
            <div>
              <div className="text-sm text-midGrey mb-1">Total saved</div>
              <div className="text-3xl font-bold text-black">{kpis.total}</div>
            </div>
            <div>
              <div className="text-sm text-midGrey mb-1">Avg stability</div>
              <div className="text-3xl font-bold text-black">{kpis.avgStability}</div>
            </div>
          </div>
        </motion.div>

        {/* Right Section - 3D Viewer (Smaller, Better Positioned) */}
        <motion.div
          className="flex justify-center items-center"
          initial={{ opacity: 0, scale: 0.95 }}
          animate={{ opacity: 1, scale: 1 }}
          transition={{ duration: 0.8, delay: 0.2 }}
        >
          <div className="w-full max-w-[450px] h-[400px] rounded-2xl overflow-hidden border border-lightGrey shadow-xl bg-white relative">
            {/* Subtle glow effect */}
            <div className="absolute inset-0 bg-gradient-to-br from-blue-50/50 to-transparent pointer-events-none" />
            
            {/* 3D Viewer */}
            <div className="relative w-full h-full z-10">
              <BenzeneGLBViewer
                mode="hero"
                height={400}
                className="w-full h-full"
              />
            </div>
            
            {/* Label */}
            <div className="absolute bottom-4 left-1/2 transform -translate-x-1/2 bg-black/80 text-white px-4 py-2 rounded-full text-sm backdrop-blur-sm">
              Molecule Preview: Benzene (C₆H₆)
            </div>
          </div>
        </motion.div>
      </div>
      
      {/* Features Section */}
      <section className="max-w-7xl mx-auto px-6">
        <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
          <Card className="p-6 hover:shadow-neon-hover transition-all duration-300">
            <div className="text-3xl mb-4">🧬</div>
            <h3 className="text-xl font-semibold text-black mb-2">AI Molecule Generator</h3>
            <p className="text-darkGrey text-sm leading-relaxed">
              Build molecules instantly using machine learning. Generate novel structures from natural language prompts.
            </p>
          </Card>
          <Card className="p-6 hover:shadow-neon-hover transition-all duration-300">
            <div className="text-3xl mb-4">🔬</div>
            <h3 className="text-xl font-semibold text-black mb-2">3D Molecular Viewer</h3>
            <p className="text-darkGrey text-sm leading-relaxed">
              Interact with high-quality 3D visualizations. Rotate, zoom, and explore molecular structures in real-time.
            </p>
          </Card>
          <Card className="p-6 hover:shadow-neon-hover transition-all duration-300">
            <div className="text-3xl mb-4">📚</div>
            <h3 className="text-xl font-semibold text-black mb-2">Your Molecule Library</h3>
            <p className="text-darkGrey text-sm leading-relaxed">
              Save, manage, and explore your creations. Build a personal collection of molecular structures.
            </p>
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
                className="min-w-[200px] flex-shrink-0"
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ duration: 0.3 }}
              >
                <div 
                  className="cursor-pointer h-full"
                  onClick={() => navigate(`/lab?id=${m.id}`)}
                >
                <Card className="overflow-hidden hover:shadow-neon-hover transition-all h-full">
                  <div className="h-32 flex items-center justify-center overflow-hidden bg-offwhite">
                    {m.molfile ? (
                      <BarbellViewer
                        molfile={m.molfile}
                        mode="card"
                        height={128}
                        atomScale={0.2}
                        bondRadius={0.05}
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
