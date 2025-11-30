import React, { useEffect } from 'react'
import { Canvas, useThree } from '@react-three/fiber'
import { OrbitControls } from '@react-three/drei'
import gsap from 'gsap'
import MoleculeScene from './MoleculeScene'
import { useKeyboardShortcuts } from './hooks/useKeyboardShortcuts'
import GridHelper from './viewer/GridHelper'
import AtomGhostMesh from './viewer/AtomGhostMesh'

/**
 * Animates camera on mount
 * Must be used inside Canvas component
 */
function CameraAnimator() {
  try {
    const { camera } = useThree()

    useEffect(() => {
      // Store initial position
      const startPos = { x: 10, y: 10, z: 10 }
      camera.position.set(startPos.x, startPos.y, startPos.z)

      // Animate camera in
      gsap.to(camera.position, {
        x: 6,
        y: 6,
        z: 6,
        duration: 1.2,
        ease: 'power3.out',
      })
    }, [camera])

    return null
  } catch (error) {
    console.error('CameraAnimator must be used inside Canvas:', error)
    return null
  }
}

function KeyboardShortcuts() {
  useKeyboardShortcuts()
  return null
}

export default function SceneRoot(){
  return (
    <div style={{width:'100%', height:'100%', background: '#ffffff'}}>
      <KeyboardShortcuts />
      <Canvas 
        camera={{ position: [10, 10, 10], fov: 45 }}
        shadows
        style={{ background: '#ffffff' }}
        gl={{ 
          antialias: true,
          alpha: false,
          preserveDrawingBuffer: false,
          powerPreference: 'high-performance',
          stencil: false,
          depth: true
        }}
        dpr={[1, 2]}
        onCreated={({ gl, scene, camera }) => {
          gl.setPixelRatio(Math.min(window.devicePixelRatio, 2))
          
          // Handle WebGL context loss - prevent default to allow recovery
          const handleContextLost = (e: Event) => {
            e.preventDefault()
            console.warn('⚠️ WebGL context lost - attempting recovery')
          }
          
          const handleContextRestored = () => {
            console.log('✅ WebGL context restored - reinitializing')
            // Force re-render by updating pixel ratio
            gl.setPixelRatio(Math.min(window.devicePixelRatio, 2))
            // Clear any stale state
            gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
          }
          
          gl.domElement.addEventListener('webglcontextlost', handleContextLost)
          gl.domElement.addEventListener('webglcontextrestored', handleContextRestored)
          
          // Cleanup on unmount
          return () => {
            gl.domElement.removeEventListener('webglcontextlost', handleContextLost)
            gl.domElement.removeEventListener('webglcontextrestored', handleContextRestored)
            // Dispose of renderer resources
            gl.dispose()
          }
        }}
      >
        <color attach="background" args={[0xffffff]} />
        <ambientLight intensity={0.7} />
        <directionalLight position={[10, 10, 5]} intensity={0.8} castShadow />
        <GridHelper />
        <React.Suspense fallback={null}>
          <MoleculeScene />
        </React.Suspense>
        <AtomGhostMesh />
        <CameraAnimator />
        <OrbitControls 
          ref={(ref) => {
            if (ref) {
              // Store ref for camera focus hook
              ;(window as any).__orbitControls = ref
            }
          }}
          makeDefault 
          enablePan={true}
          enableDamping
          dampingFactor={0.05}
          minDistance={2}
          maxDistance={50}
          rotateSpeed={0.5}
        />
      </Canvas>
    </div>
  )
}

