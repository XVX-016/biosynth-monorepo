import React from 'react'
import { motion } from 'framer-motion'
import { FiMousePointer, FiPlus, FiLink, FiTrash2 } from 'react-icons/fi'
import { useLabStore } from '../../store/labStore'
import type { ToolName } from '../../types/molecule'
import ElementPalette from './ElementPalette'
import Inspector from './Inspector'
import ExploreModeUI from './ExploreModeUI'

const tools: Array<{ id: ToolName; icon: React.ReactNode; label: string }> = [
  { id: 'select', icon: <FiMousePointer size={20} />, label: 'Select' },
  { id: 'add_atom', icon: <FiPlus size={20} />, label: 'Add Atom' },
  { id: 'bond', icon: <FiLink size={20} />, label: 'Bond' },
  { id: 'delete', icon: <FiTrash2 size={20} />, label: 'Delete' },
]

export default function LabToolPanel() {
  const currentTool = useLabStore(s => s.currentTool)
  const setTool = useLabStore(s => s.setTool)

  return (
    <motion.div
      initial={{ opacity: 0, y: 10 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.4, ease: 'easeOut' }}
      className="flex flex-col items-center py-4 gap-3 w-full"
    >
      {tools.map((t, index) => {
        const active = currentTool === t.id

        return (
          <motion.button
            key={t.id}
            initial={{ opacity: 0, x: -10 }}
            animate={{ opacity: 1, x: 0 }}
            transition={{ duration: 0.3, delay: index * 0.05 }}
            whileHover={{ scale: 1.07 }}
            whileTap={{ scale: 0.95 }}
            onClick={(e) => {
              // Ripple effect
              const circle = document.createElement('span')
              circle.className = 'ripple'
              const rect = e.currentTarget.getBoundingClientRect()
              const size = Math.max(rect.width, rect.height)
              circle.style.width = circle.style.height = `${size}px`
              e.currentTarget.appendChild(circle)
              setTimeout(() => circle.remove(), 400)
              
              setTool(t.id)
            }}
            className={`
              relative w-10 h-10 rounded-lg flex items-center justify-center border
              transition-colors overflow-hidden
              ${active
                ? 'border-[#4676ff] text-[#4676ff] bg-white shadow-sm'
                : 'border-[#ccc] text-gray-600'}
              hover:bg-white
            `}
            title={t.label}
          >
            <motion.div
              animate={{
                opacity: active ? 1 : 0.7,
                y: active ? -1 : 0,
              }}
              transition={{ duration: 0.25 }}
            >
              {t.icon}
            </motion.div>
          </motion.button>
        )
      })}

      {/* Divider */}
      <div className="w-full border-t border-[#e5e7eb] my-2" />

      {/* Element Palette - shown when add_atom tool is active */}
      {currentTool === 'add_atom' && (
        <motion.div
          initial={{ opacity: 0, height: 0 }}
          animate={{ opacity: 1, height: 'auto' }}
          exit={{ opacity: 0, height: 0 }}
          transition={{ duration: 0.3 }}
          className="w-full px-2"
        >
          <div className="text-xs text-gray-500 mb-2 text-center font-medium">Elements</div>
          <ElementPalette />
        </motion.div>
      )}

      {/* Inspector - always visible */}
      <div className="w-full px-2 mt-auto">
        <Inspector />
      </div>

      {/* Explore Mode */}
      <div className="w-full px-2 mt-3">
        <div className="text-xs text-gray-500 mb-2 text-center font-medium">Explore</div>
        <ExploreModeUI />
      </div>
    </motion.div>
  )
}
