/**
 * Toolbar - Top toolbar with tool buttons and undo/redo
 */

import React from 'react';
import { useToolMode } from '../hooks/useToolMode';
import { useLab } from '../hooks/useLab';

const tools = [
  { id: 'select', label: 'Select', icon: '👆' },
  { id: 'atom', label: 'Atom', icon: '⚛️' },
  { id: 'bond', label: 'Bond', icon: '🔗' },
  { id: 'drag', label: 'Drag', icon: '✋' },
  { id: 'eraser', label: 'Erase', icon: '🗑️' },
];

export default function Toolbar() {
  const { activeTool, setTool, is } = useToolMode();
  const { undo, redo, actionManager } = useLab();

  return (
    <div className="h-12 border-b border-neutral-700 bg-neutral-800 flex items-center px-3 gap-3">
      {tools.map((tool) => (
        <button
          key={tool.id}
          onClick={() => setTool(tool.id)}
          className={`px-3 py-1 rounded text-sm transition-colors ${
            is(tool.id)
              ? 'bg-blue-600 text-white'
              : 'bg-neutral-700 text-neutral-300 hover:bg-neutral-600'
          }`}
          title={tool.label}
        >
          <span className="mr-1">{tool.icon}</span>
          {tool.label}
        </button>
      ))}

      <div className="flex-1" />

      <button
        onClick={undo}
        disabled={!actionManager.canUndo()}
        className="px-3 py-1 rounded text-sm bg-neutral-700 text-neutral-300 hover:bg-neutral-600 disabled:opacity-50 disabled:cursor-not-allowed"
      >
        Undo
      </button>
      <button
        onClick={redo}
        disabled={!actionManager.canRedo()}
        className="px-3 py-1 rounded text-sm bg-neutral-700 text-neutral-300 hover:bg-neutral-600 disabled:opacity-50 disabled:cursor-not-allowed"
      >
        Redo
      </button>
    </div>
  );
}

