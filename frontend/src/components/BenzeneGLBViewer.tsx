/**
 * BenzeneGLBViewer - Renders a 3D GLB model of benzene
 * 
 * Loads the benzene.glb file from public/chemistry-benzene/source/
 * Supports autorotate and interactive controls
 */
import React, { useRef, Suspense } from 'react';
import { Canvas, useFrame } from '@react-three/fiber';
import { OrbitControls, useGLTF, Environment } from '@react-three/drei';
import * as THREE from 'three';

// Simple error boundary component
class ErrorBoundary extends React.Component<
  { children: React.ReactNode },
  { hasError: boolean }
> {
  constructor(props: any) {
    super(props);
    this.state = { hasError: false };
  }

  static getDerivedStateFromError() {
    return { hasError: true };
  }

  componentDidCatch(error: unknown, info: unknown) {
    console.error('BenzeneGLBViewer error:', error, info);
  }

  render() {
    if (this.state.hasError) {
      return (
        <mesh>
          <boxGeometry args={[1, 1, 1]} />
          <meshStandardMaterial color="#ff0000" />
        </mesh>
      );
    }
    return this.props.children;
  }
}

interface BenzeneGLBViewerProps {
  mode?: 'hero' | 'card';
  className?: string;
  height?: number;
}

// Auto-rotate the model
function BenzeneModel({ autorotate = false }: { autorotate?: boolean }) {
  const gltf = useGLTF('/chemistry-benzene/source/Benzene.glb');
  const groupRef = useRef<THREE.Group>(null);

  // Clone the scene to avoid issues with multiple instances
  const clonedScene = React.useMemo(() => {
    if (!gltf || !gltf.scene) {
      console.error('Failed to load GLB model');
      return null;
    }
    console.log('Benzene GLB model loaded successfully');
    return gltf.scene.clone();
  }, [gltf]);

  // Center and scale the model
  React.useEffect(() => {
    if (clonedScene) {
      const box = new THREE.Box3().setFromObject(clonedScene);
      const center = box.getCenter(new THREE.Vector3());
      const size = box.getSize(new THREE.Vector3());
      const maxDim = Math.max(size.x, size.y, size.z);
      const scale = 2 / maxDim; // Scale to fit in a 2-unit box

      clonedScene.scale.multiplyScalar(scale);
      clonedScene.position.sub(center.multiplyScalar(scale));
    }
  }, [clonedScene]);

  // Autorotate for hero mode
  useFrame((_, delta) => {
    if (autorotate && groupRef.current) {
      groupRef.current.rotation.y += delta * 0.5;
    }
  });

  if (!clonedScene) {
    return null;
  }

  return (
    <group ref={groupRef}>
      <primitive object={clonedScene} />
    </group>
  );
}

// Preload the model for better performance
useGLTF.preload('/chemistry-benzene/source/Benzene.glb');

export default function BenzeneGLBViewer({
  mode = 'hero',
  className = 'w-full h-full',
  height = 500,
}: BenzeneGLBViewerProps) {
  const autorotate = mode === 'hero';
  const interactive = mode === 'hero';

  return (
    <div
      className={`relative ${className} rounded-md overflow-hidden bg-transparent`}
      style={{ height }}
    >
      <Canvas
        camera={{
          position: [0, 0, 5],
          fov: 50,
        }}
        style={{ width: '100%', height: '100%' }}
      >
        <ambientLight intensity={0.6} />
        <directionalLight position={[5, 5, 5]} intensity={0.8} />
        <directionalLight position={[-5, -3, -2]} intensity={0.4} />
        <Environment preset="city" environmentIntensity={0.6} />
        
        <ErrorBoundary>
          <Suspense fallback={
            <mesh>
              <boxGeometry args={[1, 1, 1]} />
              <meshStandardMaterial color="#888888" />
            </mesh>
          }>
            <BenzeneModel autorotate={autorotate} />
          </Suspense>
        </ErrorBoundary>
        
        {interactive && (
          <OrbitControls
            enablePan={true}
            enableZoom={true}
            enableRotate={true}
            enableDamping={true}
            dampingFactor={0.05}
            autoRotate={autorotate}
            autoRotateSpeed={0.5}
          />
        )}
      </Canvas>
    </div>
  );
}

