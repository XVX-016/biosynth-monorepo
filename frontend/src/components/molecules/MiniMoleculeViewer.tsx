import React, { Suspense } from 'react';
import { Canvas } from '@react-three/fiber';
import { OrbitControls, Environment } from '@react-three/drei';
import * as THREE from 'three';
import { MoleculeGraph } from '@biosynth/engine';
import { moleculeToRenderable } from '../../lib/engineAdapter';

interface MiniMoleculeViewerProps {
  molecule: MoleculeGraph | null;
  className?: string;
}

const ELEMENT_RADII: Record<string, number> = {
  H: 0.3,
  C: 0.5,
  O: 0.45,
  N: 0.48,
  F: 0.45,
  S: 0.6,
  P: 0.58,
  Cl: 0.55,
  Br: 0.63,
  I: 0.7,
};

function MoleculeScene({ molecule }: { molecule: MoleculeGraph }) {
  const renderable = moleculeToRenderable(molecule);

  return (
    <>
      <ambientLight intensity={0.5} />
      <pointLight position={[5, 5, 5]} intensity={1} />
      <pointLight position={[-5, -5, -5]} intensity={0.8} />

      {renderable.atoms.map((atom) => {
        const radius = ELEMENT_RADII[atom.element] || 0.5;
        return (
          <mesh key={atom.id} position={atom.position}>
            <sphereGeometry args={[radius, 32, 32]} />
            <meshPhysicalMaterial
              color="#C0C5D2"
              metalness={0.9}
              roughness={0.2}
              envMapIntensity={1.2}
            />
          </mesh>
        );
      })}

      {renderable.bonds.map((bond) => {
        const vFrom = new THREE.Vector3(...bond.from);
        const vTo = new THREE.Vector3(...bond.to);
        const diff = new THREE.Vector3().subVectors(vTo, vFrom);
        const length = diff.length();
        const mid = new THREE.Vector3().addVectors(vFrom, vTo).multiplyScalar(0.5);
        const q = new THREE.Quaternion().setFromUnitVectors(
          new THREE.Vector3(0, 1, 0),
          diff.clone().normalize()
        );
        const radius = bond.order === 1 ? 0.1 : bond.order === 2 ? 0.14 : 0.18;

        return (
          <mesh
            key={bond.id}
            position={[mid.x, mid.y, mid.z]}
            quaternion={[q.x, q.y, q.z, q.w]}
          >
            <cylinderGeometry args={[radius, radius, length, 24]} />
            <meshPhysicalMaterial
              color="#C0C5D2"
              metalness={0.9}
              roughness={0.2}
              envMapIntensity={1.2}
            />
          </mesh>
        );
      })}

      <OrbitControls
        enableZoom={false}
        enablePan={false}
        autoRotate
        autoRotateSpeed={1}
        minPolarAngle={0.1}
        maxPolarAngle={Math.PI - 0.1}
      />

      <Environment preset="city" environmentIntensity={0.6} />
    </>
  );
}

export default function MiniMoleculeViewer({ molecule, className = '' }: MiniMoleculeViewerProps) {
  if (!molecule) {
    return (
      <div className={`flex items-center justify-center bg-offwhite rounded-lg ${className}`}>
        <div className="text-midGrey text-sm">No molecule</div>
      </div>
    );
  }

  return (
    <div className={`bg-offwhite rounded-lg overflow-hidden ${className}`}>
      <Canvas
        camera={{ position: [0, 0, 6], fov: 50 }}
        gl={{ antialias: true, alpha: true, powerPreference: 'high-performance' }}
        style={{ background: 'transparent' }}
        onCreated={({ gl }) => {
          const handleContextLost = (e: Event) => {
            e.preventDefault()
            console.warn('⚠️ WebGL context lost in MiniMoleculeViewer')
          }
          const handleContextRestored = () => {
            console.log('✅ WebGL context restored in MiniMoleculeViewer')
            gl.setPixelRatio(Math.min(window.devicePixelRatio, 2))
          }
          gl.domElement.addEventListener('webglcontextlost', handleContextLost)
          gl.domElement.addEventListener('webglcontextrestored', handleContextRestored)
          return () => {
            gl.domElement.removeEventListener('webglcontextlost', handleContextLost)
            gl.domElement.removeEventListener('webglcontextrestored', handleContextRestored)
          }
        }}
      >
        <Suspense fallback={null}>
          <MoleculeScene molecule={molecule} />
        </Suspense>
      </Canvas>
    </div>
  );
}

