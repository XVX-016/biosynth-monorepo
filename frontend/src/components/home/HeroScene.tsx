import React, { useRef, Suspense } from 'react';
import { Canvas, useFrame } from '@react-three/fiber';
import { OrbitControls, Environment, EffectComposer, ChromaticAberration } from '@react-three/drei';
import * as THREE from 'three';
import { MoleculeGraph } from '@biosynth/engine';
import AtomMesh from '../r3f/AtomMesh';
import BondMesh from '../r3f/BondMesh';
import { moleculeToRenderable } from '../../lib/engineAdapter';
import { listMolecules } from '../../lib/api';
import { createMoleculeFromTemplate, getTemplateByName } from '../../data/molecule-templates';

interface HeroSceneProps {
  molecule?: MoleculeGraph | null;
}

function FloatingMolecule({ molecule }: { molecule: MoleculeGraph }) {
  const groupRef = useRef<THREE.Group>(null);
  const renderable = moleculeToRenderable(molecule);
  
  // Floating animation
  useFrame((state) => {
    if (groupRef.current) {
      // Subtle rotation
      groupRef.current.rotation.y += 0.002;
      
      // Floating motion
      const time = state.clock.elapsedTime;
      groupRef.current.position.y = Math.sin(time * 0.5) * 0.2;
    }
  });
  
  return (
    <group ref={groupRef}>
      {/* Atoms */}
      {renderable.atoms.map((atom) => (
        <AtomMesh
          key={atom.id}
          id={atom.id}
          position={atom.position}
          element={atom.element as any}
        />
      ))}
      
      {/* Bonds */}
      {renderable.bonds.map((bond) => (
        <BondMesh
          key={bond.id}
          from={bond.from}
          to={bond.to}
          order={bond.order}
        />
      ))}
    </group>
  );
}

function SceneContent({ molecule }: { molecule: MoleculeGraph | null }) {
  if (!molecule) return null;
  
  return (
    <>
      {/* Cold white/ivory lighting - NO blue */}
      <ambientLight intensity={0.4} color="#F6F7F8" />
      <directionalLight position={[5, 5, 5]} intensity={1.2} color="#FFFFFF" />
      <directionalLight position={[-5, -5, -5]} intensity={0.6} color="#F6F7F8" />
      <pointLight position={[0, 10, 0]} intensity={0.5} color="#FFFFFF" />
      
      <FloatingMolecule molecule={molecule} />
      
      <OrbitControls
        enableZoom={false}
        enablePan={false}
        autoRotate
        autoRotateSpeed={0.5}
        minPolarAngle={Math.PI / 3}
        maxPolarAngle={Math.PI / 2.2}
      />
      
      {/* Chrome-like environment with cold white reflections */}
      <Environment
        environmentIntensity={0.8}
        environmentRotation={[0, 0, 0]}
        preset="city"
      />
      
      <EffectComposer>
        <ChromaticAberration
          offset={[0.001, 0.001]}
          radialModulation={true}
          modulationOffset={0.15}
        />
      </EffectComposer>
    </>
  );
}

export default function HeroScene({ molecule: propMolecule }: HeroSceneProps) {
  const [featuredMolecule, setFeaturedMolecule] = React.useState<MoleculeGraph | null>(
    propMolecule || null
  );
  
  React.useEffect(() => {
    // Load featured molecule
    if (propMolecule) {
      setFeaturedMolecule(propMolecule);
      return;
    }
    
    // Try to load from API
    let cancelled = false;
    (async () => {
      try {
        const molecules = await listMolecules(1);
        if (!cancelled && molecules.length > 0) {
          // Try to load from saved molecule
          // For now, use a template as fallback
          const template = getTemplateByName('Benzene');
          if (template) {
            setFeaturedMolecule(createMoleculeFromTemplate(template));
          }
        } else {
          // Use default template
          const template = getTemplateByName('Benzene') || getTemplateByName('Water');
          if (template) {
            setFeaturedMolecule(createMoleculeFromTemplate(template));
          }
        }
      } catch (error) {
        console.error('Failed to load featured molecule:', error);
        // Fallback to template
        const template = getTemplateByName('Benzene') || getTemplateByName('Water');
        if (template) {
          setFeaturedMolecule(createMoleculeFromTemplate(template));
        }
      }
    })();
    
    return () => {
      cancelled = true;
    };
  }, [propMolecule]);
  
  if (!featuredMolecule) {
    return (
      <div className="w-full h-full flex items-center justify-center bg-ionBlack">
        <div className="text-chrome">Loading molecule...</div>
      </div>
    );
  }
  
  return (
    <div className="w-full h-full bg-ionBlack relative overflow-hidden">
      {/* Subtle neonCyan edge bloom */}
      <div className="absolute inset-0 bg-gradient-to-br from-neonCyan/5 via-transparent to-plasmaTeal/5 pointer-events-none" />
      <Canvas
        camera={{ position: [0, 0, 8], fov: 50 }}
        gl={{ antialias: true, alpha: true, toneMappingExposure: 1.2 }}
        style={{ background: 'transparent' }}
      >
        <Suspense fallback={null}>
          <SceneContent molecule={featuredMolecule} />
        </Suspense>
      </Canvas>
    </div>
  );
}

