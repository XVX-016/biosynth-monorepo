import React, { useRef, useMemo, Suspense } from 'react';
import { Canvas, useFrame, useThree } from '@react-three/fiber';
import { OrbitControls, Environment } from '@react-three/drei';
import { EffectComposer, ChromaticAberration, Bloom } from '@react-three/postprocessing';
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

const ELEMENT_RADII: Record<string, number> = {
  H: 0.55,
  C: 1.0,
  O: 0.9,
  N: 0.95,
  F: 0.9,
  S: 1.2,
  P: 1.15,
  Cl: 1.1,
  Br: 1.25,
  I: 1.4,
};

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
      {/* Atoms with material override */}
      {renderable.atoms.map((atom) => {
        const radius = ELEMENT_RADII[atom.element] || 1.0;
        return (
          <mesh key={atom.id} position={atom.position}>
            <sphereGeometry args={[radius, 64, 64]} />
            <meshPhysicalMaterial
              color="#C0C5D2"
              metalness={0.9}
              roughness={0.2}
              envMapIntensity={1.5}
            />
          </mesh>
        );
      })}
      
      {/* Bonds with material override */}
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
        const radius = bond.order === 1 ? 0.14 : bond.order === 2 ? 0.18 : 0.22;
        
        return (
          <mesh
            key={bond.id}
            position={[mid.x, mid.y, mid.z]}
            quaternion={[q.x, q.y, q.z, q.w]}
          >
            <cylinderGeometry args={[radius, radius, length, 48]} />
            <meshPhysicalMaterial
              color="#C0C5D2"
              metalness={0.9}
              roughness={0.2}
              envMapIntensity={1.5}
            />
          </mesh>
        );
      })}
    </group>
  );
}

function SceneContent({ molecule }: { molecule: MoleculeGraph | null }) {
  if (!molecule) return null;
  
  return (
    <>
      {/* Ivory & Chrome lighting - NO blue */}
      <ambientLight intensity={0.4} color="#F6F7F8" />
      <pointLight position={[5, 5, 5]} intensity={1.4} color="#C0C5D2" />
      <pointLight position={[-5, -5, -5]} intensity={1.4} color="#C0C5D2" />
      
      {/* Rim-light with plasmaTeal tone */}
      <pointLight position={[0, 0, -8]} intensity={0.8} color="#3BC7C9" />
      <pointLight position={[0, 8, 0]} intensity={0.6} color="#3BC7C9" />
      
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
      
      <PostProcessingEffects />
    </>
  );
}

function PostProcessingEffects() {
  const { gl } = useThree();

  const supportsComposer = useMemo(() => {
    if (!gl) return false;
    return gl.capabilities.isWebGL2 || gl.extensions.has('WEBGL_draw_buffers');
  }, [gl]);

  const passes = useMemo(
    () => [
      <Bloom
        key="bloom"
        intensity={1.5}
        luminanceThreshold={0.9}
        luminanceSmoothing={0.9}
        height={300}
        opacity={0.8}
      />,
      <ChromaticAberration
        key="chromatic"
        offset={[0.001, 0.001]}
        radialModulation
        modulationOffset={0.15}
      />,
    ],
    []
  );

  if (!supportsComposer) {
    return null;
  }

  return <EffectComposer>{passes}</EffectComposer>;
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
    <div className="w-full h-full relative overflow-hidden frosted-glass border border-chrome/30 rounded-xl" style={{
      background: 'linear-gradient(135deg, rgba(227, 230, 235, 0.1), rgba(192, 197, 210, 0.05))',
      borderImage: 'linear-gradient(135deg, #E3E6EB, #C0C5D2) 1'
    }}>
      {/* neonCyan glow bloom */}
      <div className="absolute inset-0 bg-gradient-to-br from-neonCyan/10 via-transparent to-neonCyan/5 pointer-events-none" style={{
        boxShadow: 'inset 0 0 60px rgba(139, 243, 255, 0.2)'
      }} />
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

