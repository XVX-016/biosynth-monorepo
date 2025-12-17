import React, { Suspense, useRef } from 'react';
import { Canvas, useFrame, useThree } from '@react-three/fiber';
import { OrbitControls, Environment } from '@react-three/drei';
import { EffectComposer, ChromaticAberration, Bloom } from '@react-three/postprocessing';
import * as THREE from 'three';
import type { MoleculeGraph } from '@biosynth/engine';
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

const ELEMENT_COLORS: Record<string, string> = {
  H: '#FFFFFF',
  C: '#2C2C2C', // Dark premium grey
  O: '#FF4D4D', // Premium red
  N: '#4D79FF', // Premium blue
  Cl: '#4DFF88',
  S: '#FFD700',
};

function FloatingMolecule({ molecule }: { molecule: MoleculeGraph }) {
  const groupRef = useRef<THREE.Group>(null);
  const renderable = moleculeToRenderable(molecule);

  useFrame((state) => {
    if (!groupRef.current) return;
    groupRef.current.rotation.y += 0.002;
    const time = state.clock.elapsedTime;
    groupRef.current.position.y = Math.sin(time * 0.5) * 0.2;
  });

  return (
    <group ref={groupRef}>
      {renderable.atoms.map((atom) => {
        const radius = ELEMENT_RADII[atom.element] || 1.0;
        const color = ELEMENT_COLORS[atom.element] || '#C0C5D2';
        return (
          <mesh key={atom.id} position={atom.position}>
            <sphereGeometry args={[radius, 64, 64]} />
            <meshPhysicalMaterial
              color={color}
              metalness={atom.element === 'C' ? 0.9 : 0.4}
              roughness={atom.element === 'H' ? 0.1 : 0.3}
              transmission={atom.element === 'H' ? 0.6 : 0}
              thickness={atom.element === 'H' ? 1 : 0}
              envMapIntensity={1.5}
            />
          </mesh>
        );
      })}

      {renderable.bonds.map((bond) => {
        // compute geometry transforms per bond (kept simple and deterministic)
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
              color="#E5E5E5"
              metalness={0.8}
              roughness={0.1}
              envMapIntensity={1.2}
              transmission={0.2}
              thickness={0.5}
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
      <ambientLight intensity={0.4} color="#F6F7F8" />
      <pointLight position={[5, 5, 5]} intensity={1.4} color="#C0C5D2" />
      <pointLight position={[-5, -5, -5]} intensity={1.4} color="#C0C5D2" />
      <pointLight position={[0, 0, -8]} intensity={0.8} color="#3BC7C9" />
      <pointLight position={[0, 8, 0]} intensity={0.6} color="#3BC7C9" />

      <FloatingMolecule molecule={molecule} />

      <OrbitControls
        enableZoom={false}
        enablePan={false}
        autoRotate
        autoRotateSpeed={0.4}
        minPolarAngle={0.05}
        maxPolarAngle={Math.PI - 0.05}
        rotateSpeed={0.8}
      />

      <Environment environmentIntensity={0.8} environmentRotation={[0, 0, 0]} preset="city" />

      <PostProcessingEffects />
    </>
  );
}

// Lightweight error boundary specifically to isolate Canvas rendering errors
class CanvasErrorBoundary extends React.Component<
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
    // eslint-disable-next-line no-console
    console.error('[CanvasErrorBoundary] caught error:', error, info);
  }

  render() {
    if (this.state.hasError) {
      // render a simple fallback that is friendly and non-blocking
      return (
        <div className="w-full h-full flex items-center justify-center bg-ionBlack">
          <div className="text-chrome">Graphics unavailable — try reloading the page.</div>
        </div>
      );
    }

    return <>{this.props.children}</>;
  }
}

function PostProcessingEffects() {
  const { gl, size } = useThree();
  const [isReady, setIsReady] = React.useState(false);
  const [hasError, setHasError] = React.useState(false);

  const chromaticOffset = React.useMemo(() => new THREE.Vector2(0.0006, 0.0006), []);

  const supportsPostProcessing = React.useMemo(() => {
    try {
      if (!gl) return false;
      const isWebGL2 = !!(gl as any).capabilities?.isWebGL2;
      const extGetter = (gl as any).getExtension ? (name: string) => (gl as any).getExtension(name) : () => null;
      const hasDrawBuffers = !!extGetter('WEBGL_draw_buffers') || !!extGetter('EXT_draw_buffers');
      const validSize = !!size && size.width > 0 && size.height > 0;
      return (isWebGL2 || hasDrawBuffers) && validSize;
    } catch (e) {
      // eslint-disable-next-line no-console
      console.warn('supportsPostProcessing check failed', e);
      return false;
    }
  }, [gl, size]);

  React.useEffect(() => {
    let mounted = true;
    const t = setTimeout(() => {
      if (!mounted) return;
      setIsReady(Boolean(gl) && supportsPostProcessing);
      setHasError(!supportsPostProcessing);
    }, 50);
    return () => {
      mounted = false;
      clearTimeout(t);
    };
  }, [gl, supportsPostProcessing]);

  if (!isReady || hasError) {
    return null;
  }

  try {
    return (
      <EffectComposer multisampling={4}>
        <Bloom intensity={0.6} luminanceThreshold={0.98} luminanceSmoothing={0.6} height={200} opacity={0.45} />
        <ChromaticAberration offset={chromaticOffset} radialModulation={false} modulationOffset={0.06} />
      </EffectComposer>
    );
  } catch (err) {
    // eslint-disable-next-line no-console
    console.error('[PostProcessingEffects] render error:', err);
    setHasError(true);
    return null;
  }
}

export default function HeroScene({ molecule: propMolecule }: HeroSceneProps) {
  const [featuredMolecule, setFeaturedMolecule] = React.useState<MoleculeGraph | null>(
    propMolecule || null
  );

  React.useEffect(() => {
    if (propMolecule) {
      setFeaturedMolecule(propMolecule);
      return;
    }

    let cancelled = false;
    (async () => {
      try {
        // Try to get Caffeine first as it's the hero favorite
        const template = getTemplateByName('Caffeine');
        if (template && !cancelled) {
          setFeaturedMolecule(createMoleculeFromTemplate(template));
        } else if (!cancelled) {
          // Fallback to library fetch if templates fail
          const molecules = await listMolecules(1);
          if (molecules.length > 0) {
            const fallback = getTemplateByName('Caffeine');
            if (fallback) setFeaturedMolecule(createMoleculeFromTemplate(fallback));
          }
        }
      } catch (error) {
        console.error('Failed to load hero molecule:', error);
        const fallback = getTemplateByName('Caffeine') || getTemplateByName('Water');
        if (fallback && !cancelled) {
          setFeaturedMolecule(createMoleculeFromTemplate(fallback));
        }
      }
    })();

    return () => {
      cancelled = true;
    };
  }, [propMolecule]);

  if (!featuredMolecule) {
    return (
      <div className="w-full h-full flex items-center justify-center">
        <div className="w-8 h-8 border-2 border-indigo-600 border-t-transparent rounded-full animate-spin"></div>
      </div>
    );
  }

  return (
    <div className="w-full h-full relative overflow-hidden">
      <CanvasErrorBoundary>
        <Canvas
          camera={{ position: [0, 0, 8], fov: 50 }}
          gl={{ antialias: true, alpha: true, toneMappingExposure: 1.2 }}
          style={{ background: 'transparent' }}
        >
          <Suspense fallback={null}>
            <SceneContent molecule={featuredMolecule} />
          </Suspense>
        </Canvas>
      </CanvasErrorBoundary>
    </div>
  );
}
