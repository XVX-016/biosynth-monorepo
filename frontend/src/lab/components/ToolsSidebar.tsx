/**
 * ToolsSidebar - Left sidebar with tool buttons
 * 
 * Groups tools into logical sections:
 * - Basic Editing (cursor, atom, bond, erase, select)
 * - Chemistry Tools (periodic table, validation, cleanup)
 * - AI/Modeling (predict, generate, convert)
 */

import React from 'react'
import { toolController, ToolMode } from '../engines/ToolController'
import { commandStack } from '../engines/CommandStack'

interface ToolsSidebarProps {
  selectedElement?: string | null
  onElementSelect?: (element: string | null) => void
  onValidate?: () => void
  onPredict?: () => void
}

export default function ToolsSidebar({
  selectedElement,
  onElementSelect,
  onValidate,
  onPredict,
}: ToolsSidebarProps) {
  const activeTool = toolController.activeTool
  const canUndo = commandStack.canUndo()
  const canRedo = commandStack.canRedo()

  const ToolButton = ({ 
    mode, 
    label, 
    icon 
  }: { 
    mode: ToolMode
    label: string
    icon?: string
  }) => (
    <button
      onClick={() => toolController.setTool(mode)}
      className={`w-full p-3 text-left border-b border-gray-200 transition-colors ${
        activeTool === mode
          ? 'bg-blue-50 text-blue-700 font-semibold'
          : 'bg-white hover:bg-gray-50 text-gray-700'
      }`}
    >
      <div className="flex items-center gap-2">
        {icon && <span>{icon}</span>}
        <span>{label}</span>
      </div>
    </button>
  )

  return (
    <div className="w-64 h-full bg-white border-r border-gray-300 flex flex-col">
      {/* Header */}
      <div className="p-4 border-b border-gray-300 bg-gray-50">
        <h2 className="text-lg font-semibold text-gray-800">Tools</h2>
        <p className="text-xs text-gray-500 mt-1">Select a tool to edit</p>
      </div>

      {/* Basic Editing Tools */}
      <div className="flex-1 overflow-y-auto">
        <div className="p-2">
          <h3 className="text-xs font-semibold text-gray-500 uppercase tracking-wide mb-2 px-2">
            Basic Editing
          </h3>
          <div className="space-y-0">
            <ToolButton mode={ToolMode.cursor} label="Cursor" icon="↖" />
            <ToolButton mode={ToolMode.select} label="Select" icon="✓" />
            <ToolButton mode={ToolMode.addAtom} label="Add Atom" icon="⚛" />
            <ToolButton mode={ToolMode.addBond} label="Add Bond" icon="─" />
            <ToolButton mode={ToolMode.erase} label="Erase" icon="✕" />
          </div>
        </div>

        {/* Chemistry Tools */}
        <div className="p-2 border-t border-gray-200">
          <h3 className="text-xs font-semibold text-gray-500 uppercase tracking-wide mb-2 px-2">
            Chemistry
          </h3>
          <div className="space-y-0">
            <button
              onClick={() => {
                // TODO: Open periodic table modal
                console.log('Periodic table')
              }}
              className="w-full p-3 text-left border-b border-gray-200 bg-white hover:bg-gray-50 text-gray-700"
            >
              Periodic Table
            </button>
            <button
              onClick={onValidate}
              className="w-full p-3 text-left border-b border-gray-200 bg-white hover:bg-gray-50 text-gray-700"
            >
              Validate Structure
            </button>
            <button
              onClick={() => {
                // TODO: Cleanup geometry
                console.log('Cleanup geometry')
              }}
              className="w-full p-3 text-left border-b border-gray-200 bg-white hover:bg-gray-50 text-gray-700"
            >
              Cleanup Geometry
            </button>
          </div>
        </div>

        {/* AI/Modeling */}
        <div className="p-2 border-t border-gray-200">
          <h3 className="text-xs font-semibold text-gray-500 uppercase tracking-wide mb-2 px-2">
            AI & Modeling
          </h3>
          <div className="space-y-0">
            <button
              onClick={onPredict}
              className="w-full p-3 text-left border-b border-gray-200 bg-white hover:bg-gray-50 text-gray-700"
            >
              Predict Properties
            </button>
            <button
              onClick={() => {
                // TODO: Generate molecule
                console.log('Generate molecule')
              }}
              className="w-full p-3 text-left border-b border-gray-200 bg-white hover:bg-gray-50 text-gray-700"
            >
              Generate Molecule
            </button>
            <button
              onClick={() => {
                // TODO: Convert format
                console.log('Convert format')
              }}
              className="w-full p-3 text-left border-b border-gray-200 bg-white hover:bg-gray-50 text-gray-700"
            >
              Convert Format
            </button>
          </div>
        </div>
      </div>

      {/* Undo/Redo */}
      <div className="p-2 border-t border-gray-300 bg-gray-50">
        <div className="flex gap-2">
          <button
            onClick={() => commandStack.undo()}
            disabled={!canUndo}
            className={`flex-1 p-2 text-sm rounded ${
              canUndo
                ? 'bg-white hover:bg-gray-100 text-gray-700 border border-gray-300'
                : 'bg-gray-100 text-gray-400 cursor-not-allowed'
            }`}
          >
            Undo
          </button>
          <button
            onClick={() => commandStack.redo()}
            disabled={!canRedo}
            className={`flex-1 p-2 text-sm rounded ${
              canRedo
                ? 'bg-white hover:bg-gray-100 text-gray-700 border border-gray-300'
                : 'bg-gray-100 text-gray-400 cursor-not-allowed'
            }`}
          >
            Redo
          </button>
        </div>
      </div>
    </div>
  )
}

