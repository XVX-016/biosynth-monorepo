import React from 'react'
import { motion, AnimatePresence } from 'framer-motion'
import { Element } from '@biosynth/engine'
import { useMoleculeStore } from '../store/moleculeStore'

const ELEMENTS: Array<{ element: Element; label: string; color: string }> = [
  { element: 'H', label: 'H', color: '#3b82f6' },
  { element: 'C', label: 'C', color: '#1e293b' },
  { element: 'N', label: 'N', color: '#8b5cf6' },
  { element: 'O', label: 'O', color: '#ef4444' },
  { element: 'F', label: 'F', color: '#10b981' },
  { element: 'S', label: 'S', color: '#f59e0b' },
  { element: 'P', label: 'P', color: '#f97316' },
  { element: 'Cl', label: 'Cl', color: '#84cc16' },
  { element: 'Br', label: 'Br', color: '#991b1b' },
  { element: 'I', label: 'I', color: '#7c3aed' },
]

export default function AtomPalette() {
  const tool = useMoleculeStore((state) => state.tool)
  const atomToAdd = useMoleculeStore((state) => state.atomToAdd)
  const setAtomToAdd = useMoleculeStore((state) => state.setAtomToAdd)

  if (tool !== 'add-atom') return null

  return (
    <AnimatePresence>
      <motion.div
        initial={{ opacity: 0, y: -20 }}
        animate={{ opacity: 1, y: 0 }}
        exit={{ opacity: 0, y: -20 }}
        className="fixed left-20 top-4 z-50 frosted-glass rounded-lg shadow-glass border border-chrome/20 p-3 backdrop-blur-md"
      >
        <div className="text-xs text-chrome mb-2">Select Element</div>
        <div className="grid grid-cols-5 gap-2">
          {ELEMENTS.map((el) => (
            <motion.button
              key={el.element}
              whileHover={{ scale: 1.1 }}
              whileTap={{ scale: 0.95 }}
              onClick={() => setAtomToAdd(el.element)}
              className={`w-10 h-10 rounded-lg font-semibold text-white transition-all ${
                atomToAdd === el.element
                  ? 'ring-2 ring-neonCyan ring-offset-2 shadow-neon-sm'
                  : ''
              }`}
              style={{ backgroundColor: el.color }}
              title={el.label}
            >
              {el.label}
            </motion.button>
          ))}
        </div>
      </motion.div>
    </AnimatePresence>
  )
}

