/**
 * LabLayout - Complete Lab Page Layout
 * 
 * Phase 10: Lab Page UI Rebuild (Complete)
 * 
 * Layout:
 * - Left: Toolbar + Element Palette
 * - Center: Editor (2D/3D tabs)
 * - Right: Prediction Panel + Inspector
 * - Bottom: Console (validation messages)
 */

import React, { useState } from 'react'
import { EditorProvider } from './EditorContext'
import { EditorPanel } from './panels/EditorPanel'
import { ToolbarPanel } from './panels/ToolbarPanel'
import { InspectorPanel } from './panels/InspectorPanel'
import { ConsolePanel } from './panels/ConsolePanel'
import { PredictionPanel } from './PredictionPanel'
import { Molecule } from '@/lib/molecule'
import { useEditorContext } from './EditorContext'

interface LabLayoutProps {
  initialMolecule?: Molecule
}

function LabLayoutContent() {
  const { molecule } = useEditorContext()
  const [rightPanelCollapsed, setRightPanelCollapsed] = useState(false)

  return (
      <div className="h-screen flex flex-col bg-gray-50">
        {/* Main content area */}
        <div className="flex-1 flex overflow-hidden">
          {/* Left sidebar */}
          <div className="w-64 bg-white border-r border-gray-200 p-4 overflow-y-auto">
            <ToolbarPanel />
            
            {/* Element palette placeholder */}
            <div className="mt-4 bg-white border border-gray-200 rounded-lg p-3">
              <h3 className="text-xs font-semibold text-gray-700 mb-2">Elements</h3>
              <div className="grid grid-cols-5 gap-2">
                {['C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'H'].map((el) => (
                  <button
                    key={el}
                    className="px-2 py-1.5 text-xs rounded border border-gray-300 bg-white text-gray-700 hover:bg-gray-50"
                  >
                    {el}
                  </button>
                ))}
              </div>
            </div>
          </div>

          {/* Center: Editor */}
          <div className="flex-1 flex flex-col overflow-hidden">
            <div className="flex-1 p-4">
              <EditorPanel />
            </div>
            
            {/* Bottom: Console */}
            <div className="p-4 border-t border-gray-200 bg-white">
              <ConsolePanel />
            </div>
          </div>

          {/* Right sidebar */}
          <div className={`w-80 bg-white border-l border-gray-200 transition-all duration-300 ${
            rightPanelCollapsed ? 'w-0 overflow-hidden' : 'p-4 overflow-y-auto'
          }`}>
            {!rightPanelCollapsed && (
              <div className="space-y-4">
                <div className="flex items-center justify-between">
                  <h2 className="text-sm font-semibold text-gray-700">Properties</h2>
                  <button
                    onClick={() => setRightPanelCollapsed(true)}
                    className="text-xs text-gray-500 hover:text-gray-700"
                  >
                    ←
                  </button>
                </div>
                
                <PredictionPanel 
                  molecule={useEditorContext().molecule}
                  className=""
                />
                
                <InspectorPanel />
              </div>
            )}
          </div>

          {/* Collapsed right panel button */}
          {rightPanelCollapsed && (
            <button
              onClick={() => setRightPanelCollapsed(false)}
              className="absolute right-0 top-1/2 transform -translate-y-1/2 bg-white border border-gray-200 rounded-l-lg px-2 py-4 text-xs text-gray-500 hover:text-gray-700"
            >
              →
            </button>
          )}
        </div>
      </div>
  )
}

export function LabLayout({ initialMolecule }: LabLayoutProps) {
  return (
    <EditorProvider initialMolecule={initialMolecule}>
      <LabLayoutContent />
    </EditorProvider>
  )
}


