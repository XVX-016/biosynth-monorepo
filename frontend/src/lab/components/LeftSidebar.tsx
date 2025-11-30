/**
 * LeftSidebar - Left sidebar with tools and inspector
 */

import React from 'react';
import ToolsPanel from './panels/ToolsPanel';
import InspectorPanel from './panels/InspectorPanel';
import ElementPalette from './ElementPalette';
import ValidationPanel from './ValidationPanel';

export default function LeftSidebar() {
  return (
    <div className="w-64 h-full bg-neutral-800 border-r border-neutral-700 flex flex-col overflow-y-auto">
      <ElementPalette />
      <ToolsPanel />
      <div className="flex-1">
        <InspectorPanel />
      </div>
      <ValidationPanel />
    </div>
  );
}

