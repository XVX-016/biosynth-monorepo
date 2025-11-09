import React from 'react'
import { motion } from 'framer-motion'
import { useMoleculeStore, Tool } from '../store/moleculeStore'
import { useHistoryStore } from '../store/historyStore'
import { undo, redo } from '../store/historyStore'

export default function ToolPanel() {
  const tool = useMoleculeStore((state) => state.tool)
  const setTool = useMoleculeStore((state) => state.setTool)
  const currentBondOrder = useMoleculeStore((state) => state.currentBondOrder)
  const setBondOrder = useMoleculeStore((state) => state.setBondOrder)
  const canUndo = useHistoryStore((state) => state.canUndo)
  const canRedo = useHistoryStore((state) => state.canRedo)

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
      className="fixed left-4 top-1/2 -translate-y-1/2 z-50 bg-panel rounded-lg shadow-elev-1 p-3 space-y-2"
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
              ? 'bg-accent-blue text-white'
              : 'bg-aluminum-DEFAULT text-text-primary hover:bg-aluminum-dark'
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
          <div className="text-xs text-text-secondary text-center mb-1">Bond Order</div>
          {[1, 2, 3].map((order) => (
            <motion.button
              key={order}
              whileHover={{ scale: 1.05 }}
              whileTap={{ scale: 0.95 }}
              onClick={() => setBondOrder(order)}
              className={`w-full px-2 py-1 text-xs rounded ${
                currentBondOrder === order
                  ? 'bg-accent-blue text-white'
                  : 'bg-aluminum-DEFAULT text-text-primary'
              }`}
            >
              {order === 1 ? 'Single' : order === 2 ? 'Double' : 'Triple'}
            </motion.button>
          ))}
        </motion.div>
      )}

      {/* Undo/Redo buttons */}
      <div className="pt-2 border-t border-aluminum-dark space-y-1">
        <motion.button
          whileHover={{ scale: 1.05 }}
          whileTap={{ scale: 0.95 }}
          onClick={undo}
          disabled={!canUndo}
          className={`w-12 h-10 rounded-lg text-sm ${
            canUndo
              ? 'bg-aluminum-DEFAULT text-text-primary hover:bg-aluminum-dark'
              : 'bg-aluminum-light text-text-tertiary cursor-not-allowed'
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
          className={`w-12 h-10 rounded-lg text-sm ${
            canRedo
              ? 'bg-aluminum-DEFAULT text-text-primary hover:bg-aluminum-dark'
              : 'bg-aluminum-light text-text-tertiary cursor-not-allowed'
          }`}
          title="Redo (Ctrl+Y)"
        >
          ↷
        </motion.button>
      </div>
    </motion.div>
  )
}

