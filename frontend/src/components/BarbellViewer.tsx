/**
 * BarbellViewer - Ball-and-stick 3D molecule viewer using react-three-fiber
 * 
 * Supports two modes:
 * - "card": Lazy mount, no autorotate, minimal controls (for library cards)
 * - "hero": Always mounted, autorotate, interactive (for homepage)
 */
import React, { useRef, useMemo, useState } from 'react';
import { Canvas, useFrame } from '@react-three/fiber';
import { OrbitControls } from '@react-three/drei';
import * as THREE from 'three';
import { parseMolfile } from '../utils/molfileParser';
import { getAtomColorNumeric, DEFAULT_ATOM_COLOR_NUMERIC } from '../utils/atomColors';

interface BarbellViewerProps {
  molfile?: string | null;
  mode?: 'card' | 'hero' | 'preview';
  className?: string;
  atomScale?: number;
  bondRadius?: number;
  height?: number;
  interactive?: boolean; // Enable OrbitControls when true
  autorotate?: boolean; // Auto-rotate when true
  hovered?: boolean; // Hover state for card mode (enables interaction)
  highlightAtoms?: number[]; // Atom indices to highlight (for KAB, spectroscopy, etc.)
  highlightColor?: number; // Color for highlighted atoms (default: yellow)
}

// Use shared color system
const ELEMENT_COLORS: Record<string, number> = {
  // Core supported elements use shared colors
  H: getAtomColorNumeric('H'),
  C: getAtomColorNumeric('C'),
  N: getAtomColorNumeric('N'),
  O: getAtomColorNumeric('O'),
  F: getAtomColorNumeric('F'),
  S: getAtomColorNumeric('S'),
  Cl: getAtomColorNumeric('Cl'),
  Br: getAtomColorNumeric('Br'),
  // Extended elements (for compatibility with existing molecules)
  P: 0xff8000,   // Orange
  I: 0x940094,   // Purple
  Li: 0xcc80ff,  // Violet
  Na: 0xab5cf2,  // Violet
  K: 0x8f40d4,   // Violet
  Mg: 0x8aff00,  // Green
  Ca: 0x3dff00,  // Green
  Fe: 0xe06633,  // Orange
  Cu: 0xc88033,  // Orange
  Zn: 0x7d80b0,  // Blue-gray
};

const DEFAULT_COLOR = DEFAULT_ATOM_COLOR_NUMERIC;

// Bond cylinder component
function BondCylinder({
  start,
  end,
  radius,
  quality = 'high'
}: {
  start: THREE.Vector3;
  end: THREE.Vector3;
  radius: number;
  quality?: 'high' | 'low';
}) {
  const length = start.distanceTo(end);
  const mid = new THREE.Vector3().addVectors(start, end).multiplyScalar(0.5);

  // Calculate rotation to align cylinder
  const direction = new THREE.Vector3().subVectors(end, start).normalize();
  const yAxis = new THREE.Vector3(0, 1, 0);
  const quaternion = new THREE.Quaternion().setFromUnitVectors(yAxis, direction);

  // Reduced segments for low quality (card mode) to improve performance
  const segments = quality === 'low' ? 8 : 16;

  return (
    <mesh position={mid} quaternion={quaternion}>
      <cylinderGeometry args={[radius, radius, length, segments]} />
      <meshStandardMaterial color={0xcccccc} metalness={0.25} roughness={0.5} />
    </mesh>
  );
}

// Main molecule scene component
function MoleculeScene({
  molfile,
  autorotate = false,
  isUserInteracting = false,
  atomScale = 0.25,
  bondRadius = 0.06,
  quality = 'high',
  highlightAtoms = [],
  highlightColor = 0xffff00,
}: {
  molfile: string;
  autorotate?: boolean;
  isUserInteracting?: boolean;
  atomScale?: number;
  bondRadius?: number;
  quality?: 'high' | 'low';
  highlightAtoms?: number[];
  highlightColor?: number;
}) {
  const groupRef = useRef<THREE.Group>(null);

  const { atoms, bonds, centroid } = useMemo(() => {
    if (!molfile) {
      return { atoms: [], bonds: [], centroid: new THREE.Vector3(0, 0, 0) };
    }

    try {
      const parsed = parseMolfile(molfile);
      const atoms = parsed.atoms;
      const bonds = parsed.bonds;

      console.log('BarbellViewer parsed:', { atomCount: atoms.length, bondCount: bonds.length });

      // Calculate centroid
      if (atoms.length === 0) {
        console.warn('BarbellViewer: No atoms parsed from molfile');
        return { atoms: [], bonds: [], centroid: new THREE.Vector3(0, 0, 0) };
      }

      const centroid = atoms.reduce(
        (acc, atom) => {
          acc.x += atom.x;
          acc.y += atom.y;
          acc.z += atom.z;
          return acc;
        },
        { x: 0, y: 0, z: 0 }
      );
      centroid.x /= atoms.length;
      centroid.y /= atoms.length;
      centroid.z /= atoms.length;

      return {
        atoms,
        bonds,
        centroid: new THREE.Vector3(centroid.x, centroid.y, centroid.z),
      };
    } catch (e) {
      console.warn('Molfile parse failed', e);
      return { atoms: [], bonds: [], centroid: new THREE.Vector3(0, 0, 0) };
    }
  }, [molfile]);

  // Autorotate for hero mode - only Y-axis (sideways) and pause when user is interacting
  useFrame((_, delta) => {
    if (autorotate && !isUserInteracting && groupRef.current) {
      // Only rotate on Y-axis (sideways rotation)
      groupRef.current.rotation.y += delta * 0.5;
    }
  });

  if (atoms.length === 0) {
    return (
      <group ref={groupRef}>
        <mesh rotation={[0.6, 0.2, 0]}>
          <icosahedronGeometry args={[0.6, 0]} />
          <meshStandardMaterial color={0x7aa6ff} metalness={0.2} roughness={0.3} />
        </mesh>
      </group>
    );
  }

  // Prepare bond vectors (centered)
  const bondVectors = useMemo(() => {
    return bonds
      .map((bond) => {
        const a1 = atoms[bond.a1];
        const a2 = atoms[bond.a2];
        if (!a1 || !a2) return null;

        const start = new THREE.Vector3(a1.x, a1.y, a1.z).sub(centroid);
        const end = new THREE.Vector3(a2.x, a2.y, a2.z).sub(centroid);
        return { start, end };
      })
      .filter((v): v is { start: THREE.Vector3; end: THREE.Vector3 } => v !== null);
  }, [atoms, bonds, centroid]);

  return (
    <group ref={groupRef} position={centroid.clone().multiplyScalar(-1)}>
      {/* Render atoms */}
      {atoms.map((atom, i) => {
        const pos = new THREE.Vector3(atom.x, atom.y, atom.z).sub(centroid);
        // Reduced segments for low quality (card mode) to improve performance
        const segments = quality === 'low' ? 12 : 32;
        const isHighlighted = highlightAtoms.includes(i);
        const baseColor = ELEMENT_COLORS[atom.element] || DEFAULT_COLOR;
        const color = isHighlighted ? highlightColor : baseColor;
        const finalScale = isHighlighted ? atomScale * 1.3 : atomScale;
        return (
          <mesh key={`atom-${i}`} position={pos.toArray()}>
            <sphereGeometry args={[finalScale, segments, segments]} />
            <meshStandardMaterial
              color={color}
              metalness={isHighlighted ? 0.5 : 0.2}
              roughness={isHighlighted ? 0.3 : 0.4}
              emissive={isHighlighted ? highlightColor : 0x000000}
              emissiveIntensity={isHighlighted ? 0.3 : 0}
            />
          </mesh>
        );
      })}

      {/* Render bonds */}
      {bondVectors.map((bond, i) => (
        <BondCylinder
          key={`bond-${i}`}
          start={bond.start}
          end={bond.end}
          radius={bondRadius}
          quality={quality}
        />
      ))}
    </group>
  );
}

// Error boundary for Canvas rendering errors
class CanvasErrorBoundary extends React.Component<
  { children: React.ReactNode; fallback?: React.ReactNode },
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
    console.warn('[BarbellViewer] Canvas error:', error, info);
  }

  render() {
    if (this.state.hasError) {
      return this.props.fallback || (
        <div className="w-full h-full flex items-center justify-center bg-gray-50 text-gray-400 text-xs">
          Preview unavailable
        </div>
      );
    }
    return <>{this.props.children}</>;
  }
}

// Main component
export default function BarbellViewer({
  molfile,
  mode = 'card',
  className = 'w-full h-full',
  atomScale = 0.25,
  bondRadius = 0.06,
  height = 200,
  interactive: interactiveProp,
  autorotate: autorotateProp,
  hovered = false,
  highlightAtoms = [],
  highlightColor = 0xffff00,
}: BarbellViewerProps) {
  const [isUserInteracting, setIsUserInteracting] = useState(false);
  const [contextLost, setContextLost] = useState(false);
  const rendererRef = useRef<THREE.WebGLRenderer | null>(null);

  const autorotate = autorotateProp !== undefined ? autorotateProp : (mode === 'hero');
  // For card mode, only enable interaction when hovered
  const interactive = interactiveProp !== undefined
    ? interactiveProp
    : (mode === 'hero' || (mode === 'card' && hovered));
  // Use high quality for all modes - better visual quality
  const quality = 'high';

  // CRITICAL: Cleanup WebGL resources on unmount
  React.useEffect(() => {
    return () => {
      if (rendererRef.current) {
        try {
          // Dispose of all resources
          rendererRef.current.dispose();
          // Force WebGL context loss to free GPU memory
          rendererRef.current.forceContextLoss();
          rendererRef.current = null;
        } catch (e) {
          console.warn('[BarbellViewer] Cleanup error:', e);
        }
      }
    };
  }, []);

  // Calculate camera distance based on molecule size
  const cameraRadius = useMemo(() => {
    if (!molfile) return 4;

    try {
      const { atoms } = parseMolfile(molfile);
      if (!atoms.length) return 4;

      // Calculate bounding sphere
      const center = atoms.reduce(
        (acc, a) => {
          acc.x += a.x;
          acc.y += a.y;
          acc.z += a.z;
          return acc;
        },
        { x: 0, y: 0, z: 0 }
      );
      center.x /= atoms.length;
      center.y /= atoms.length;
      center.z /= atoms.length;

      let maxDist = 0;
      atoms.forEach((a) => {
        const d = Math.hypot(a.x - center.x, a.y - center.y, a.z - center.z);
        if (d > maxDist) maxDist = d;
      });

      return Math.max(2.5, maxDist * 3.0);
    } catch {
      return 4;
    }
  }, [molfile]);

  // Fallback when no molfile
  if (!molfile) {
    return (
      <div
        className={`relative ${className} bg-gray-50 flex items-center justify-center overflow-hidden rounded-md`}
        style={{ height }}
      >
        <div className="text-gray-400 text-sm">No molecule data</div>
      </div>
    );
  }

  // Fallback when WebGL context is lost
  if (contextLost) {
    return (
      <div
        className={`relative ${className} bg-gray-50 flex items-center justify-center overflow-hidden rounded-md`}
        style={{ height }}
      >
        <div className="text-gray-400 text-xs text-center px-2">
          Preview unavailable<br />
          <span className="text-[10px] opacity-60">WebGL context lost</span>
        </div>
      </div>
    );
  }

  return (
    <CanvasErrorBoundary
      fallback={
        <div
          className={`relative ${className} bg-gray-50 flex items-center justify-center overflow-hidden rounded-md`}
          style={{ height }}
        >
          <div className="text-gray-400 text-xs">Preview unavailable</div>
        </div>
      }
    >
      <div
        className={`relative ${className} rounded-md overflow-hidden`}
        style={{ height }}
      >
        <Canvas
          camera={{
            position: [0, 0, Math.max(4, cameraRadius)],
            fov: 45,
          }}
          style={{ width: '100%', height: '100%', background: '#ffffff' }}
          dpr={Math.min(window.devicePixelRatio, 2)}
          performance={{ min: 0.5 }}
          flat={false}
          shadows={mode !== 'card'}
          gl={{
            antialias: mode !== 'card',
            alpha: true,
            preserveDrawingBuffer: false,
            powerPreference: 'high-performance',
            stencil: false,
            depth: true,
          }}
          onCreated={({ gl }) => {
            // Store renderer reference for cleanup
            rendererRef.current = gl;

            // Handle WebGL context loss
            gl.domElement.addEventListener('webglcontextlost', (e) => {
              e.preventDefault();
              console.warn('[BarbellViewer] WebGL context lost');
              setContextLost(true);
            });
            gl.domElement.addEventListener('webglcontextrestored', () => {
              console.log('[BarbellViewer] WebGL context restored');
              setContextLost(false);
            });
          }}
        >

          <ambientLight intensity={0.6} />
          <directionalLight position={[5, 5, 5]} intensity={0.9} />
          <directionalLight position={[-5, -3, -2]} intensity={0.5} />
          <directionalLight position={[0, 5, -5]} intensity={0.3} />
          <MoleculeScene
            molfile={molfile}
            autorotate={autorotate}
            isUserInteracting={isUserInteracting}
            atomScale={atomScale}
            bondRadius={bondRadius}
            quality={quality}
            highlightAtoms={highlightAtoms}
            highlightColor={highlightColor}
          />
          {interactive && (
            <OrbitControls
              enablePan={false}
              enableZoom={false}
              enableRotate={true}
              enableDamping={true}
              dampingFactor={0.05}
              autoRotate={false}
              minDistance={mode === 'card' ? 2 : 1}
              maxDistance={mode === 'card' ? 8 : 20}
              onStart={() => setIsUserInteracting(true)}
              onEnd={() => setIsUserInteracting(false)}
            />
          )}
        </Canvas>
      </div>
    </CanvasErrorBoundary>
  );
}
