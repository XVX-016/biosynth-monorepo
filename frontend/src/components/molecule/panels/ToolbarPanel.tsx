/**
 * ToolbarPanel - Tool selection toolbar
 * 
 * Phase 10: Lab Page UI Rebuild
 * Phase 11: Will be expanded with full toolbar
 */

import React from 'react'
import { useEditorContext } from '../EditorContext'
import type { EditorTool } from '@/lib/molecule'

const tools: { id: EditorTool; label: string; icon: string }[] = [
  { id: 'select', label: 'Select', icon: '👆' },
  { id: 'add-atom', label: 'Atom', icon: '⚛️' },
  { id: 'bond', label: 'Bond', icon: '🔗' },
  { id: 'delete', label: 'Erase', icon: '🗑️' },
  { id: 'move', label: 'Move', icon: '↔️' },
]

export function ToolbarPanel() {
  const { tool, setTool, canUndo, canRedo, undo, redo } = useEditorContext()

  return (
    <div className="bg-white border border-gray-200 rounded-lg p-3">
      <h3 className="text-xs font-semibold text-gray-700 mb-3">Tools</h3>
      
      <div className="grid grid-cols-3 gap-2 mb-4">
        {tools.map((t) => (
          <button
            key={t.id}
            onClick={() => setTool(t.id)}
            className={`px-3 py-2 text-xs rounded border ${
              tool === t.id
                ? 'bg-blue-50 border-blue-500 text-blue-700'
                : 'bg-white border-gray-300 text-gray-700 hover:bg-gray-50'
            }`}
            title={t.label}
          >
            <div className="text-lg mb-1">{t.icon}</div>
            <div>{t.label}</div>
          </button>
        ))}
      </div>

      <div className="flex gap-2">
        <button
          onClick={undo}
          disabled={!canUndo}
          className="flex-1 px-3 py-1.5 text-xs rounded border border-gray-300 bg-white text-gray-700 disabled:opacity-50 disabled:cursor-not-allowed hover:bg-gray-50"
        >
          ↶ Undo
        </button>
        <button
          onClick={redo}
          disabled={!canRedo}
          className="flex-1 px-3 py-1.5 text-xs rounded border border-gray-300 bg-white text-gray-700 disabled:opacity-50 disabled:cursor-not-allowed hover:bg-gray-50"
        >
          ↷ Redo
        </button>
      </div>
    </div>
  )
}

