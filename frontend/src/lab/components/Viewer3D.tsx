/**
 * Viewer3D - 3D molecule viewer using 3Dmol.js
 * 
 * Displays the molecule in 3D and stays in sync with 2D editor
 */

import { useEffect, useRef } from 'react'
import { moleculeEngine } from '../engines/MoleculeStateEngine'

interface Viewer3DProps {
  width?: number
  height?: number
  style?: 'stick' | 'ballstick' | 'sphere' | 'cartoon'
}

declare global {
  interface Window {
    $3Dmol: any
  }
}

export default function Viewer3D({ 
  width = 800, 
  height = 600,
  style = 'ballstick'
}: Viewer3DProps) {
  const containerRef = useRef<HTMLDivElement>(null)
  const viewerRef = useRef<any>(null)

  // Load 3Dmol.js if not already loaded
  useEffect(() => {
    if (window.$3Dmol) {
      initViewer()
      return
    }

    // Load 3Dmol.js from CDN
    const script = document.createElement('script')
    script.src = 'https://3Dmol.org/build/3Dmol-min.js'
    script.async = true
    script.onload = () => {
      initViewer()
    }
    document.body.appendChild(script)

    return () => {
      // Cleanup
      if (viewerRef.current) {
        viewerRef.current.clear()
      }
    }
  }, [])

  // Initialize 3Dmol viewer
  const initViewer = () => {
    if (!containerRef.current || !window.$3Dmol) return

    // Create viewer
    viewerRef.current = window.$3Dmol.createViewer(containerRef.current, {
      defaultcolors: window.$3Dmol.rasmolElementColors,
    })

    // Set background
    viewerRef.current.setBackgroundColor(0xffffff)

    // Load initial molecule
    updateMolecule()
  }

  // Update molecule when engine changes
  const updateMolecule = () => {
    if (!viewerRef.current || !window.$3Dmol) return

    try {
      const sdf = moleculeEngine.toSDF()
      
      if (!sdf || sdf.trim().length === 0) {
        viewerRef.current.clear()
        return
      }

      // Clear previous model
      viewerRef.current.clear()

      // Add model
      viewerRef.current.addModel(sdf, 'sdf')

      // Set style based on prop
      const styleConfig: Record<string, any> = {
        stick: { stick: {} },
        ballstick: { stick: {}, sphere: { scale: 0.2 } },
        sphere: { sphere: { scale: 0.3 } },
        cartoon: { cartoon: {} },
      }

      viewerRef.current.setStyle({}, styleConfig[style] || styleConfig.ballstick)

      // Zoom to fit
      viewerRef.current.zoomTo()

      // Render
      viewerRef.current.render()
    } catch (error) {
      console.error('Error updating 3D viewer:', error)
    }
  }

  // Watch for molecule changes
  useEffect(() => {
    if (!viewerRef.current) return

    // Update when molecule changes
    // This is a simple polling approach - in production, you'd use a proper state management
    const interval = setInterval(() => {
      updateMolecule()
    }, 500) // Update every 500ms

    return () => clearInterval(interval)
  }, [style])

  // Update when style changes
  useEffect(() => {
    if (viewerRef.current) {
      updateMolecule()
    }
  }, [style])

  return (
    <div 
      ref={containerRef} 
      className="w-full h-full bg-gray-900"
      style={{ width, height }}
    />
  )
}

