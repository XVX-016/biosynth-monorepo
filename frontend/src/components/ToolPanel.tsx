import React from 'react'
import { motion, AnimatePresence } from 'framer-motion'
import { useMoleculeStore, Tool } from '../store/moleculeStore'
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
      initial={{ opacity: 0, x: -20 }}
      animate={{ opacity: 1, x: 0 }}
      className="fixed left-4 top-1/2 -translate-y-1/2 z-50 frosted-glass rounded-lg shadow-glass border border-chrome/20 p-3 space-y-2 backdrop-blur-md"
    >
      {/* Tool buttons */}
      {tools.map((t) => (
        <motion.button
          key={t.id}
          whileHover={{ scale: 1.1 }}
          whileTap={{ scale: 0.95 }}
          onClick={() => setTool(t.id)}
          className={`w-12 h-12 rounded-lg font-medium transition-all ${
            tool === t.id
              ? 'bg-plasma-neon text-ionBlack shadow-neon-sm'
              : 'bg-frostedGlass text-chrome hover:text-ivory hover:border-neonCyan/30 border border-chrome/20'
          }`}
          title={`${t.label} (${t.hotkey})`}
        >
          <span className="text-xl">{t.icon}</span>
        </motion.button>
      ))}

      {/* Bond order selector (only when bond tool is active) */}
      {tool === 'bond' && (
        <motion.div
          initial={{ opacity: 0, height: 0 }}
          animate={{ opacity: 1, height: 'auto' }}
          className="mt-2 space-y-1"
        >
          <div className="text-xs text-chrome text-center mb-1">Bond Order</div>
          {[1, 2, 3].map((order) => (
            <motion.button
              key={order}
              whileHover={{ scale: 1.05 }}
              whileTap={{ scale: 0.95 }}
              onClick={() => setBondOrder(order)}
              className={`w-full px-2 py-1 text-xs rounded transition-all ${
                currentBondOrder === order
                  ? 'bg-plasma-neon text-ionBlack shadow-neon-sm'
                  : 'bg-frostedGlass text-chrome hover:text-ivory border border-chrome/20'
              }`}
            >
              {order === 1 ? 'Single' : order === 2 ? 'Double' : 'Triple'}
            </motion.button>
          ))}
        </motion.div>
      )}

      {/* Undo/Redo buttons */}
      <div className="pt-2 border-t border-chrome/20 space-y-1">
        <motion.button
          whileHover={{ scale: 1.05 }}
          whileTap={{ scale: 0.95 }}
          onClick={undo}
          disabled={!canUndo}
          className={`w-12 h-10 rounded-lg text-sm transition-all ${
            canUndo
              ? 'bg-frostedGlass text-chrome hover:text-ivory hover:border-neonCyan/30 border border-chrome/20'
              : 'bg-frostedGlass/50 text-chrome/50 cursor-not-allowed border border-chrome/10'
          }`}
          title="Undo (Ctrl+Z)"
        >
          ↶
        </motion.button>
        <motion.button
          whileHover={{ scale: 1.05 }}
          whileTap={{ scale: 0.95 }}
          onClick={redo}
          disabled={!canRedo}
          className={`w-12 h-10 rounded-lg text-sm transition-all ${
            canRedo
              ? 'bg-frostedGlass text-chrome hover:text-ivory hover:border-neonCyan/30 border border-chrome/20'
              : 'bg-frostedGlass/50 text-chrome/50 cursor-not-allowed border border-chrome/10'
          }`}
          title="Redo (Ctrl+Y)"
        >
          ↷
        </motion.button>
      </div>

      {/* Auto-bond toggle */}
      <div className="pt-2 border-t border-chrome/20 space-y-1">
        <div className="text-xs text-chrome text-center mb-1">Auto-Bond</div>
        <button
          onClick={() => setAutoBond(!autoBond)}
          className={`w-24 mx-auto h-8 rounded-lg text-sm transition-all ${
            autoBond 
              ? 'bg-plasma-neon text-ionBlack shadow-neon-sm' 
              : 'bg-frostedGlass text-chrome border border-chrome/20'
          }`}
          aria-pressed={autoBond}
          aria-label="Toggle auto-bond mode"
        >
          {autoBond ? 'On' : 'Off'}
        </button>
      </div>

      {/* Create Bond button (when bond tool is active) */}
      {tool === 'bond' && (
        <motion.div
          initial={{ opacity: 0, height: 0 }}
          animate={{ opacity: 1, height: 'auto' }}
          className="pt-2 border-t border-chrome/20"
        >
          <motion.button
            whileHover={{ scale: 1.05 }}
            whileTap={{ scale: 0.95 }}
            onClick={handleCreateBond}
            className="w-full px-3 py-2 rounded-lg text-sm transition-all bg-frostedGlass text-chrome hover:text-ivory hover:border-neonCyan/30 border border-chrome/20"
            title="Click two atoms to create a bond"
          >
            Create Bond
          </motion.button>
          <div className="text-xs text-chrome/70 text-center mt-1">
            Select 2 atoms
          </div>
        </motion.div>
      )}

      {/* Templates Section */}
      <div className="pt-2 border-t border-chrome/20">
        <motion.button
          whileHover={{ scale: 1.02 }}
          whileTap={{ scale: 0.98 }}
          onClick={() => setTemplatesExpanded(!templatesExpanded)}
          className="w-full px-3 py-2 rounded-lg text-sm transition-all bg-frostedGlass text-chrome hover:text-ivory hover:border-neonCyan/30 border border-chrome/20 flex items-center justify-between"
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

