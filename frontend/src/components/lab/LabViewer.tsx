import React from 'react'
import SceneRoot from './SceneRoot'

/**
 * Main Lab Viewer component - wraps the 3D scene with animations
 */
export default function LabViewer() {
  return (
    <div className="w-full h-full">
      <SceneRoot />
    </div>
  )
}

