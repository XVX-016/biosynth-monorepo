import React from 'react'
import { motion, AnimatePresence } from 'framer-motion'
import { useMoleculeStore } from '../store/moleculeStore'
import type { Tool } from '../store/moleculeStore'
import { useHistoryStore } from '../store/historyStore'
import { undo, redo } from '../store/historyStore'
import { createBondSafe } from '../kernel/bonds'
import { pushState } from '../store/historyStore'
import { TemplatePanel } from './TemplatePanel/TemplatePanel'

export default function ToolPanel() {
  const tool = useMoleculeStore((state) => state.tool)
  const setTool = useMoleculeStore((state) => state.setTool)
  const currentBondOrder = useMoleculeStore((state) => state.currentBondOrder)
  const setBondOrder = useMoleculeStore((state) => state.setBondOrder)
  const autoBond = useMoleculeStore((state) => state.autoBond)
  const setAutoBond = useMoleculeStore((state) => state.setAutoBond)
  const canUndo = useHistoryStore((state) => state.canUndo)
  const canRedo = useHistoryStore((state) => state.canRedo)
  const currentMolecule = useMoleculeStore((state) => state.currentMolecule)
  const selectedAtomId = useMoleculeStore((state) => state.selectedAtomId)
  const [templatesExpanded, setTemplatesExpanded] = React.useState(false)
  
  // Handle Create Bond button - creates bond between two selected atoms
  // This is a helper button that shows instructions when bond tool is active
  const handleCreateBond = () => {
    // The actual bond creation happens via BondTool when two atoms are clicked
    // This button just provides visual feedback
    if (!currentMolecule || !selectedAtomId) {
      // Switch to bond tool if not already active
      setTool('bond')
    }
  }

  const tools: Array<{ id: Tool; label: string; icon: string; hotkey?: string }> = [
    { id: 'select', label: 'Select', icon: '👆', hotkey: 'S' },
    { id: 'add-atom', label: 'Add Atom', icon: '➕', hotkey: 'A' },
    { id: 'bond', label: 'Bond', icon: '🔗', hotkey: 'B' },
    { id: 'delete', label: 'Delete', icon: '🗑️', hotkey: 'DEL' },
  ]

  // Handle keyboard shortcuts
  React.useEffect(() => {
    const handleKeyDown = (e: KeyboardEvent) => {
      // Tool shortcuts
      if (e.key === 's' || e.key === 'S') {
        if (!e.ctrlKey && !e.metaKey) {
          e.preventDefault()
          setTool('select')
        }
      } else if (e.key === 'a' || e.key === 'A') {
        if (!e.ctrlKey && !e.metaKey) {
          e.preventDefault()
          setTool('add-atom')
        }
      } else if (e.key === 'b' || e.key === 'B') {
        if (!e.ctrlKey && !e.metaKey) {
          e.preventDefault()
          setTool('bond')
        }
      } else if (e.key === 'Delete' || e.key === 'Backspace') {
        if (!e.ctrlKey && !e.metaKey) {
          e.preventDefault()
          setTool('delete')
        }
      }

      // Undo/Redo
      if ((e.ctrlKey || e.metaKey) && e.key === 'z' && !e.shiftKey) {
        e.preventDefault()
        undo()
      } else if ((e.ctrlKey || e.metaKey) && (e.key === 'y' || (e.key === 'z' && e.shiftKey))) {
        e.preventDefault()
        redo()
      }
    }

    window.addEventListener('keydown', handleKeyDown)
    return () => window.removeEventListener('keydown', handleKeyDown)
  }, [setTool])

  return (
    <motion.div
      initial={{ opacity: 0, y: 10 }}
      animate={{ opacity: 1, y: 0 }}
      className="space-y-4"
    >
      <div className="grid grid-cols-2 sm:grid-cols-4 gap-2">
        {tools.map((t) => (
          <motion.button
            key={t.id}
            whileHover={{ scale: 1.02 }}
            whileTap={{ scale: 0.98 }}
            onClick={() => setTool(t.id)}
            className={`h-12 rounded-lg font-medium transition-all border ${
              tool === t.id
                ? 'bg-black text-white border-black shadow-neon-sm'
                : 'bg-white text-darkGrey border-lightGrey hover:border-darkGrey/40 hover:text-black'
            }`}
            title={`${t.label} (${t.hotkey})`}
          >
            <span className="text-base">{t.icon}</span>
            <span className="sr-only">{t.label}</span>
          </motion.button>
        ))}
      </div>

      {/* Bond order selector (only when bond tool is active) */}
      {tool === 'bond' && (
        <motion.div
          initial={{ opacity: 0, height: 0 }}
          animate={{ opacity: 1, height: 'auto' }}
          className="space-y-2"
        >
          <div className="text-xs uppercase tracking-[0.3em] text-midGrey">Bond order</div>
          <div className="grid grid-cols-3 gap-2">
            {[1, 2, 3].map((order) => (
              <button
                key={order}
                onClick={() => setBondOrder(order)}
                className={`rounded-lg px-2 py-2 text-xs font-semibold transition-all border ${
                  currentBondOrder === order
                    ? 'bg-black text-white border-black shadow-neon-sm'
                    : 'bg-white text-darkGrey border-lightGrey hover:border-darkGrey/40 hover:text-black'
                }`}
              >
                {order === 1 ? 'Single' : order === 2 ? 'Double' : 'Triple'}
              </button>
            ))}
          </div>
        </motion.div>
      )}

      {/* Undo/Redo buttons */}
      <div className="grid grid-cols-2 gap-2">
        <motion.button
          whileHover={{ scale: 1.02 }}
          whileTap={{ scale: 0.98 }}
          onClick={undo}
          disabled={!canUndo}
          className={`h-11 rounded-lg text-sm transition-all border ${
            canUndo
              ? 'bg-white text-darkGrey border-lightGrey hover:border-darkGrey/40 hover:text-black'
              : 'bg-offwhite text-midGrey border-lightGrey cursor-not-allowed'
          }`}
          title="Undo (Ctrl+Z)"
        >
          ↶ Undo
        </motion.button>
        <motion.button
          whileHover={{ scale: 1.02 }}
          whileTap={{ scale: 0.98 }}
          onClick={redo}
          disabled={!canRedo}
          className={`h-11 rounded-lg text-sm transition-all border ${
            canRedo
              ? 'bg-white text-darkGrey border-lightGrey hover:border-darkGrey/40 hover:text-black'
              : 'bg-offwhite text-midGrey border-lightGrey cursor-not-allowed'
          }`}
          title="Redo (Ctrl+Y)"
        >
          ↷ Redo
        </motion.button>
      </div>

      {/* Auto-bond toggle */}
      <div className="flex items-center justify-between rounded-lg border border-lightGrey px-3 py-2">
        <div>
          <p className="text-xs uppercase tracking-[0.3em] text-midGrey">Auto-bond</p>
          <p className="text-sm text-black">{autoBond ? 'Enabled' : 'Disabled'}</p>
        </div>
        <button
          onClick={() => setAutoBond(!autoBond)}
          className={`inline-flex h-8 items-center rounded-full border px-4 text-sm font-medium transition ${
            autoBond ? 'bg-black text-white border-black' : 'bg-white text-darkGrey border-lightGrey'
          }`}
        >
          {autoBond ? 'On' : 'Off'}
        </button>
      </div>

      {/* Create Bond button (when bond tool is active) */}
      {tool === 'bond' && (
        <motion.div
          initial={{ opacity: 0, height: 0 }}
          animate={{ opacity: 1, height: 'auto' }}
          className="rounded-lg border border-lightGrey px-3 py-2 bg-offwhite"
        >
          <motion.button
            whileHover={{ scale: 1.01 }}
            whileTap={{ scale: 0.99 }}
            onClick={handleCreateBond}
            className="w-full rounded-md bg-white px-3 py-2 text-sm font-medium text-darkGrey border border-lightGrey"
            title="Click two atoms to create a bond"
          >
            Create Bond
          </motion.button>
          <div className="text-xs text-midGrey text-center mt-1">
            Select 2 atoms
          </div>
        </motion.div>
      )}

      {/* Templates Section */}
      <div>
        <motion.button
          whileHover={{ scale: 1.02 }}
          whileTap={{ scale: 0.98 }}
          onClick={() => setTemplatesExpanded(!templatesExpanded)}
          className="w-full px-3 py-2 rounded-lg text-sm transition-all border border-lightGrey bg-white text-darkGrey hover:border-darkGrey/40 flex items-center justify-between"
        >
          <span>Templates</span>
          <span className="text-xs">{templatesExpanded ? '▼' : '▶'}</span>
        </motion.button>
        <AnimatePresence>
          {templatesExpanded && (
            <motion.div
              initial={{ opacity: 0, height: 0 }}
              animate={{ opacity: 1, height: 'auto' }}
              exit={{ opacity: 0, height: 0 }}
              className="mt-2 overflow-hidden"
            >
              <div className="max-h-96 overflow-y-auto">
                <TemplatePanel />
              </div>
            </motion.div>
          )}
        </AnimatePresence>
      </div>
    </motion.div>
  )
}

