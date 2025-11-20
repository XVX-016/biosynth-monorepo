/**
 * BarbellViewer - Ball-and-stick 3D molecule viewer using react-three-fiber
 * 
 * Supports two modes:
 * - "card": Lazy mount, no autorotate, minimal controls (for library cards)
 * - "hero": Always mounted, autorotate, interactive (for homepage)
 */
import React, { useRef, useMemo } from 'react';
import { Canvas, useFrame } from '@react-three/fiber';
import { OrbitControls } from '@react-three/drei';
import * as THREE from 'three';
import { parseMolfile, type Atom, type Bond } from '../utils/molfileParser';

interface BarbellViewerProps {
  molfile?: string | null;
  mode?: 'card' | 'hero';
  className?: string;
  atomScale?: number;
  bondRadius?: number;
  height?: number;
}

// CPK element colors
const ELEMENT_COLORS: Record<string, number> = {
  H: 0xffffff,   // White
  C: 0x909090,   // Dark gray
  N: 0x3050f8,   // Blue
  O: 0xff0d0d,   // Red
  F: 0x90e050,   // Green
  P: 0xff8000,   // Orange
  S: 0xffff30,   // Yellow
  Cl: 0x1ff01f,  // Green
  Br: 0xa62929,  // Dark red
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

const DEFAULT_COLOR = 0xcccccc;

// Atom sphere component
function AtomSphere({ atom, scale }: { atom: Atom; scale: number }) {
  const color = ELEMENT_COLORS[atom.element] || DEFAULT_COLOR;
  
  return (
    <mesh position={[atom.x, atom.y, atom.z]}>
      <sphereGeometry args={[scale, 16, 16]} />
      <meshStandardMaterial color={color} metalness={0.15} roughness={0.5} />
    </mesh>
  );
}

// Bond cylinder component
function BondCylinder({ 
  start, 
  end, 
  radius 
}: { 
  start: THREE.Vector3; 
  end: THREE.Vector3; 
  radius: number;
}) {
  const length = start.distanceTo(end);
  const mid = new THREE.Vector3().addVectors(start, end).multiplyScalar(0.5);
  
  // Calculate rotation to align cylinder
  const direction = new THREE.Vector3().subVectors(end, start).normalize();
  const yAxis = new THREE.Vector3(0, 1, 0);
  const quaternion = new THREE.Quaternion().setFromUnitVectors(yAxis, direction);
  
  return (
    <mesh position={mid} quaternion={quaternion}>
      <cylinderGeometry args={[radius, radius, length, 8]} />
      <meshStandardMaterial color={0xcccccc} metalness={0.2} roughness={0.6} />
    </mesh>
  );
}

// Main molecule scene component
function MoleculeScene({
  molfile,
  autorotate = false,
  atomScale = 0.25,
  bondRadius = 0.06,
}: {
  molfile: string;
  autorotate?: boolean;
  atomScale?: number;
  bondRadius?: number;
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

      // Calculate centroid
      if (atoms.length === 0) {
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

  // Autorotate for hero mode
  useFrame((_, delta) => {
    if (autorotate && groupRef.current) {
      groupRef.current.rotation.y += delta * 0.35;
      groupRef.current.rotation.x += delta * 0.06;
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
        return (
          <mesh key={`atom-${i}`} position={pos.toArray()}>
            <sphereGeometry args={[atomScale, 24, 24]} />
            <meshStandardMaterial
              color={ELEMENT_COLORS[atom.element] || DEFAULT_COLOR}
              metalness={0.15}
              roughness={0.5}
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
        />
      ))}
    </group>
  );
}

// Main component
export default function BarbellViewer({
  molfile,
  mode = 'card',
  className = 'w-full h-full',
  atomScale = 0.25,
  bondRadius = 0.06,
  height = 200,
}: BarbellViewerProps) {
  const autorotate = mode === 'hero';
  const interactive = mode === 'hero';

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

  return (
    <div
      className={`relative ${className} rounded-md overflow-hidden bg-transparent`}
      style={{ height }}
    >
      <Canvas
        camera={{
          position: [0, 0, Math.max(4, cameraRadius)],
          fov: 45,
        }}
        style={{ width: '100%', height: '100%' }}
      >
        <ambientLight intensity={0.6} />
        <directionalLight position={[5, 5, 5]} intensity={0.8} />
        <directionalLight position={[-5, -3, -2]} intensity={0.4} />
        <MoleculeScene
          molfile={molfile}
          autorotate={autorotate}
          atomScale={atomScale}
          bondRadius={bondRadius}
        />
        {interactive && (
          <OrbitControls
            enablePan={true}
            enableZoom={true}
            enableRotate={true}
            enableDamping={true}
            dampingFactor={0.05}
          />
        )}
      </Canvas>
    </div>
  );
}
