import React from 'react'
import { motion, AnimatePresence } from 'framer-motion'
import type { Element } from '@biosynth/engine'
import { useMoleculeStore } from '../store/moleculeStore'

// Element colors using Ivory & Chrome palette - distinct but theme-consistent
const ELEMENTS: Array<{ element: Element; label: string; color: string }> = [
  { element: 'H', label: 'H', color: '#8BF3FF' }, // neonCyan
  { element: 'C', label: 'C', color: '#C0C5D2' }, // chrome
  { element: 'N', label: 'N', color: '#C6BDFE' }, // violetEdge
  { element: 'O', label: 'O', color: '#3BC7C9' }, // plasmaTeal
  { element: 'F', label: 'F', color: '#8BF3FF' }, // neonCyan (lighter variant)
  { element: 'S', label: 'S', color: '#F6F7F8' }, // ivory
  { element: 'P', label: 'P', color: '#C0C5D2' }, // chrome (darker)
  { element: 'Cl', label: 'Cl', color: '#3BC7C9' }, // plasmaTeal
  { element: 'Br', label: 'Br', color: '#2B2E33' }, // spaceGrey
  { element: 'I', label: 'I', color: '#C6BDFE' }, // violetEdge
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
              className={`w-10 h-10 rounded-lg font-semibold text-ionBlack transition-all ${
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

