/**
 * LabPage - Main Lab page component
 */

import React from 'react';
import LeftSidebar from './LeftSidebar';
import Toolbar from './Toolbar';
import Editor2D from './Editor2D';
import Viewer3D from './Viewer3D';
import BottomDock from './BottomDock';
import { useLab } from '../hooks/useLab';

export default function LabPage() {
  const { currentMolecule } = useLab();

  return (
    <div className="flex h-screen w-full bg-neutral-900 text-white overflow-hidden">
      <LeftSidebar />
      
      <div className="flex flex-col flex-1 relative">
        <Toolbar />
        
        <div className="flex flex-1 overflow-hidden">
          <Editor2D />
          <Viewer3D />
        </div>
        
        <BottomDock />
      </div>
    </div>
  );
}

