/**
 * ToolsPanel - Panel showing available tools
 */

import React from 'react';
import { useToolMode } from '../../hooks/useToolMode';

export default function ToolsPanel() {
  const { activeTool, setTool, is } = useToolMode();

  return (
    <div className="p-3 border-b border-neutral-700">
      <h2 className="text-sm font-semibold mb-2 opacity-70">Tools</h2>
      <div className="flex flex-col gap-2">
        <button
          onClick={() => setTool('select')}
          className={`px-2 py-1 rounded text-left text-sm transition-colors ${
            is('select')
              ? 'bg-blue-600 text-white'
              : 'bg-neutral-700 text-neutral-300 hover:bg-neutral-600'
          }`}
        >
          Select
        </button>
        <button
          onClick={() => setTool('atom')}
          className={`px-2 py-1 rounded text-left text-sm transition-colors ${
            is('atom')
              ? 'bg-blue-600 text-white'
              : 'bg-neutral-700 text-neutral-300 hover:bg-neutral-600'
          }`}
        >
          Add Atom
        </button>
        <button
          onClick={() => setTool('bond')}
          className={`px-2 py-1 rounded text-left text-sm transition-colors ${
            is('bond')
              ? 'bg-blue-600 text-white'
              : 'bg-neutral-700 text-neutral-300 hover:bg-neutral-600'
          }`}
        >
          Add Bond
        </button>
        <button
          onClick={() => setTool('drag')}
          className={`px-2 py-1 rounded text-left text-sm transition-colors ${
            is('drag')
              ? 'bg-blue-600 text-white'
              : 'bg-neutral-700 text-neutral-300 hover:bg-neutral-600'
          }`}
        >
          Drag
        </button>
        <button
          onClick={() => setTool('eraser')}
          className={`px-2 py-1 rounded text-left text-sm transition-colors ${
            is('eraser')
              ? 'bg-blue-600 text-white'
              : 'bg-neutral-700 text-neutral-300 hover:bg-neutral-600'
          }`}
        >
          Erase
        </button>
      </div>
    </div>
  );
}

