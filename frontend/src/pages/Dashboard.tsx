import React from 'react'
import { motion } from 'framer-motion'
import MoleculeViewer from '../components/MoleculeViewer'

export default function Dashboard(){
  return (
    <motion.div
      initial={{ opacity: 0 }}
      animate={{ opacity: 1 }}
      transition={{ duration: 0.3 }}
      className="space-y-6"
    >
      <header className="flex items-center justify-between">
        <div className="text-xl font-semibold text-text-primary">BioSynth AI</div>
      </header>
      <section className="grid grid-cols-3 gap-6">
        <motion.div
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ delay: 0.1 }}
          className="col-span-2 p-6 bg-panel rounded-lg shadow-elev-1"
        >
          <h2 className="text-3xl font-bold">
            Design the Future of <span className="text-accent-blue">Molecular Science</span>
          </h2>
          <p className="text-text-secondary mt-2">
            Generate, analyze, and optimize molecular structures with cutting-edge AI.
          </p>
          <div className="mt-6">
            <MoleculeViewer/>
          </div>
        </motion.div>
        <motion.aside
          initial={{ opacity: 0, x: 20 }}
          animate={{ opacity: 1, x: 0 }}
          transition={{ delay: 0.2 }}
          className="p-6 bg-panel rounded-lg shadow-elev-1 space-y-4"
        >
          <motion.button
            whileHover={{ scale: 1.05 }}
            whileTap={{ scale: 0.95 }}
            className="w-full bg-accent-blue text-white px-4 py-3 rounded-lg font-medium"
          >
            Generate Molecule
          </motion.button>
          <motion.button
            whileHover={{ scale: 1.05 }}
            whileTap={{ scale: 0.95 }}
            className="w-full bg-aluminum-DEFAULT px-4 py-3 rounded-lg font-medium text-text-primary"
          >
            Explore Library
          </motion.button>
        </motion.aside>
      </section>
    </motion.div>
  )
}

