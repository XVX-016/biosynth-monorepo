import React, { useState } from 'react'
import { motion } from 'framer-motion'
import Navbar from '../../components/Navbar'
import LabViewer from '../../components/lab/LabViewer'
import ToolButtons from '../../components/lab/ToolButtons'
import { useLabStore } from '../../store/labStore'
import ElementPalette from '../../components/lab/ElementPalette'
import Inspector from '../../components/lab/Inspector'
import ExploreModeUI from '../../components/lab/ExploreModeUI'
import PrecisionPanel from '../../components/lab/sidebar/PrecisionPanel'
import ValidationOverlay from '../../components/lab/overlay/ValidationOverlay'
import ContextMenu from '../../components/lab/ui/ContextMenu'
import { useContextMenu } from '../../components/lab/interaction/useContextMenu'
import SpectroscopyPanel from '../../components/lab/panels/SpectroscopyPanel'

/**
 * Full-screen Lab layout:
 * - Navbar at top
 * - Top toolbar below navbar
 * - Left sidebar toolbar
 * - Full workspace area (no card limitation)
 */
export default function NewLabLayout() {
  const [mobileOpen, setMobileOpen] = useState(false)
  const currentTool = useLabStore(s => s.currentTool)
  const currentElement = useLabStore(s => s.currentElement)
  const mol = useLabStore(s => s.molecule)
  const showGrid = useLabStore(s => s.showGrid ?? true)
  const setShowGrid = useLabStore(s => s.setShowGrid)
  const showValidation = useLabStore(s => s.showValidation ?? true)
  const setShowValidation = useLabStore(s => s.setShowValidation)
  const { menuState, showMenu, hideMenu, deleteAtom, deleteBond } = useContextMenu()

  return (
    <div className="h-screen w-screen flex flex-col bg-white overflow-hidden">
      {/* Standard Navbar */}
      <Navbar onToggleMenu={() => setMobileOpen((v) => !v)} />

      {/* Top Toolbar */}
      <motion.div
        initial={{ y: -10, opacity: 0 }}
        animate={{ y: 0, opacity: 1 }}
        transition={{ duration: 0.3 }}
        className="h-12 border-b border-[#e5e7eb] bg-white flex items-center justify-between px-4 shrink-0"
      >
        <div className="flex items-center gap-4">
          <span className="text-sm font-medium text-gray-700">Molecule Lab</span>
          <div className="h-4 w-px bg-gray-300" />
          <span className="text-xs text-gray-500">
            Atoms: <strong className="text-gray-700">{mol.atoms.length}</strong> | 
            Bonds: <strong className="text-gray-700">{mol.bonds.length}</strong>
          </span>
        </div>
        <div className="flex items-center gap-2">
          <label className="flex items-center gap-1.5 text-xs cursor-pointer">
            <input
              type="checkbox"
              checked={showGrid}
              onChange={(e) => setShowGrid(e.target.checked)}
              className="rounded"
            />
            <span className="text-gray-600">Grid</span>
          </label>
          <label className="flex items-center gap-1.5 text-xs cursor-pointer">
            <input
              type="checkbox"
              checked={showValidation}
              onChange={(e) => setShowValidation(e.target.checked)}
              className="rounded"
            />
            <span className="text-gray-600">Validation</span>
          </label>
          <div className="h-4 w-px bg-gray-300" />
          <button
            className="px-3 py-1.5 text-xs rounded-lg border border-gray-300 bg-white hover:bg-gray-50 transition-colors"
            onClick={() => useLabStore.getState().undo()}
          >
            Undo
          </button>
          <button
            className="px-3 py-1.5 text-xs rounded-lg border border-gray-300 bg-white hover:bg-gray-50 transition-colors"
            onClick={() => useLabStore.getState().redo()}
          >
            Redo
          </button>
          <button
            className="px-3 py-1.5 text-xs rounded-lg border border-red-300 bg-white hover:bg-red-50 text-red-600 transition-colors"
            onClick={() => useLabStore.getState().resetMolecule()}
          >
            Reset
          </button>
        </div>
      </motion.div>

      {/* Main Content Area */}
      <div className="flex flex-1 overflow-hidden">
        {/* Left Sidebar Toolbar */}
        <motion.aside
          initial={{ x: -20, opacity: 0 }}
          animate={{ x: 0, opacity: 1 }}
          transition={{ duration: 0.4, ease: 'easeOut' }}
          className="w-20 border-r border-[#e5e7eb] bg-[#f9fafb] flex flex-col shrink-0 overflow-y-auto"
        >
          <ToolButtons />
          
          {/* Element Palette - shown when add_atom tool is active */}
          {currentTool === 'add_atom' && (
            <motion.div
              initial={{ opacity: 0, height: 0 }}
              animate={{ opacity: 1, height: 'auto' }}
              exit={{ opacity: 0, height: 0 }}
              transition={{ duration: 0.3 }}
              className="w-full px-2 mt-4"
            >
              <div className="text-xs text-gray-500 mb-2 text-center font-medium">Elements</div>
              <ElementPalette />
            </motion.div>
          )}

          {/* Precision Panel */}
          <div className="w-full px-2 mt-4">
            <PrecisionPanel />
          </div>

          {/* Spectroscopy Panel */}
          <div className="w-full px-2 mt-4">
            <SpectroscopyPanel />
          </div>

          {/* Inspector - always visible */}
          <div className="w-full px-2 mt-auto mb-4">
            <Inspector />
          </div>

          {/* Explore Mode */}
          <div className="w-full px-2 mb-4">
            <div className="text-xs text-gray-500 mb-2 text-center font-medium">Explore</div>
            <ExploreModeUI />
          </div>
        </motion.aside>

        {/* Full Workspace Area - No Card Limitation */}
        <motion.main
          initial={{ opacity: 0 }}
          animate={{ opacity: 1 }}
          transition={{ duration: 0.5 }}
          className="flex-1 bg-white relative overflow-hidden"
        >
          <LabViewer />
          <ValidationOverlay />
        </motion.main>
      </div>

      {/* Context Menu */}
      <ContextMenu
        visible={menuState.visible}
        x={menuState.x}
        y={menuState.y}
        atomId={menuState.atomId}
        bondId={menuState.bondId}
        onDeleteAtom={deleteAtom}
        onDeleteBond={deleteBond}
        onClose={hideMenu}
      />
    </div>
  )
}
