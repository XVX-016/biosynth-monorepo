/**
 * Viewer3D - 3D molecule viewer using react-three-fiber
 */

import React, { useEffect, useRef } from 'react';
import { Canvas, useFrame } from '@react-three/fiber';
import { OrbitControls, PerspectiveCamera } from '@react-three/drei';
import * as THREE from 'three';
import { useLab } from '../hooks/useLab';
import { useSelection } from '../hooks/useSelection';
import type { Atom, Bond } from '../engines/MoleculeStateEngine';

// Atom sphere component
function AtomSphere({
  atom,
  selected,
  invalid,
}: {
  atom: Atom;
  selected: boolean;
  invalid: boolean;
}) {
  const color = getElementColor(atom.element);
  const scale = selected ? 1.2 : invalid ? 1.1 : 1.0;
  const emissive = invalid ? 0x440000 : selected ? 0x000044 : 0x000000;

  return (
    <mesh position={[atom.x / 20, atom.y / 20, (atom.z || 0) / 20]}>
      <sphereGeometry args={[0.3 * scale, 16, 16]} />
      <meshStandardMaterial
        color={color}
        emissive={emissive}
        emissiveIntensity={invalid ? 0.3 : selected ? 0.2 : 0}
        metalness={0.3}
        roughness={0.4}
      />
    </mesh>
  );
}

// Bond cylinder component
function BondCylinder({
  atom1,
  atom2,
  order,
  selected,
  invalid,
}: {
  atom1: Atom;
  atom2: Atom;
  order: number;
  selected: boolean;
  invalid: boolean;
}) {
  const start = new THREE.Vector3(atom1.x / 20, atom1.y / 20, (atom1.z || 0) / 20);
  const end = new THREE.Vector3(atom2.x / 20, atom2.y / 20, (atom2.z || 0) / 20);
  const length = start.distanceTo(end);
  const mid = new THREE.Vector3().addVectors(start, end).multiplyScalar(0.5);

  const direction = new THREE.Vector3().subVectors(end, start).normalize();
  const yAxis = new THREE.Vector3(0, 1, 0);
  const quaternion = new THREE.Quaternion().setFromUnitVectors(yAxis, direction);

  const color = invalid ? 0xef4444 : selected ? 0x3b82f6 : 0x666666;
  const radius = order === 2 ? 0.08 : order === 3 ? 0.1 : 0.06;

  return (
    <mesh position={mid} quaternion={quaternion}>
      <cylinderGeometry args={[radius, radius, length, 8]} />
      <meshStandardMaterial
        color={color}
        metalness={0.25}
        roughness={0.5}
        emissive={invalid ? 0x440000 : selected ? 0x000044 : 0x000000}
        emissiveIntensity={invalid ? 0.2 : selected ? 0.1 : 0}
      />
    </mesh>
  );
}

// Element colors (CPK)
function getElementColor(element: string): number {
  const colors: Record<string, number> = {
    H: 0xffffff,
    C: 0x909090,
    N: 0x3050f8,
    O: 0xff0d0d,
    F: 0x90e050,
    P: 0xff8000,
    S: 0xffff30,
    Cl: 0x1ff01f,
    Br: 0xa62929,
    I: 0x940094,
  };
  return colors[element] || 0xcccccc;
}

// Molecule scene component
function MoleculeScene({
  atoms,
  bonds,
  selectedAtomId,
  selectedBondId,
  invalidAtomIds,
  invalidBondIds,
}: {
  atoms: Atom[];
  bonds: Bond[];
  selectedAtomId: string | null;
  selectedBondId: string | null;
  invalidAtomIds: Set<string>;
  invalidBondIds: Set<string>;
}) {
  const groupRef = useRef<THREE.Group>(null);

  // Auto-rotate
  useFrame((_, delta) => {
    if (groupRef.current) {
      groupRef.current.rotation.y += delta * 0.3;
    }
  });

  // Calculate centroid for centering
  const centroid = React.useMemo(() => {
    if (atoms.length === 0) return new THREE.Vector3(0, 0, 0);
    const sum = atoms.reduce(
      (acc, atom) => {
        acc.x += atom.x / 20;
        acc.y += atom.y / 20;
        acc.z += (atom.z || 0) / 20;
        return acc;
      },
      { x: 0, y: 0, z: 0 }
    );
    return new THREE.Vector3(
      sum.x / atoms.length,
      sum.y / atoms.length,
      sum.z / atoms.length
    );
  }, [atoms]);

  return (
    <group ref={groupRef} position={centroid.clone().multiplyScalar(-1)}>
      {/* Render bonds */}
      {bonds.map((bond) => {
        const atom1 = atoms.find((a) => a.id === bond.atoms[0]);
        const atom2 = atoms.find((a) => a.id === bond.atoms[1]);
        if (!atom1 || !atom2) return null;

        return (
          <BondCylinder
            key={bond.id}
            atom1={atom1}
            atom2={atom2}
            order={bond.order}
            selected={bond.id === selectedBondId}
            invalid={invalidBondIds.has(bond.id)}
          />
        );
      })}

      {/* Render atoms */}
      {atoms.map((atom) => (
        <AtomSphere
          key={atom.id}
          atom={atom}
          selected={atom.id === selectedAtomId}
          invalid={invalidAtomIds.has(atom.id)}
        />
      ))}
    </group>
  );
}

export default function Viewer3D() {
  const { currentMolecule, validationResult } = useLab();
  const { selectedAtomId, selectedBondId } = useSelection();

  const invalidAtomIds = React.useMemo(() => {
    if (!validationResult) return new Set<string>();
    return new Set(
      validationResult.errors
        .filter((e) => e.atomId)
        .map((e) => e.atomId!)
    );
  }, [validationResult]);

  const invalidBondIds = React.useMemo(() => {
    if (!validationResult) return new Set<string>();
    return new Set(
      validationResult.errors
        .filter((e) => e.bondId)
        .map((e) => e.bondId!)
    );
  }, [validationResult]);

  if (!currentMolecule || currentMolecule.atoms.size === 0) {
    return (
      <div className="w-[32%] bg-neutral-900 border-l border-neutral-700 relative flex items-center justify-center">
        <div className="text-neutral-500 text-sm">No molecule loaded</div>
      </div>
    );
  }

  const atoms = Array.from(currentMolecule.atoms.values());
  const bonds = Array.from(currentMolecule.bonds.values());

  return (
    <div className="w-[32%] bg-neutral-900 border-l border-neutral-700 relative">
      <Canvas
        camera={{ position: [0, 0, 10], fov: 50 }}
        gl={{ antialias: true, alpha: false }}
        style={{ background: '#1a1a1a' }}
      >
        <ambientLight intensity={0.6} />
        <directionalLight position={[5, 5, 5]} intensity={0.9} />
        <directionalLight position={[-5, -3, -2]} intensity={0.5} />
        <directionalLight position={[0, 5, -5]} intensity={0.3} />

        <MoleculeScene
          atoms={atoms}
          bonds={bonds}
          selectedAtomId={selectedAtomId}
          selectedBondId={selectedBondId}
          invalidAtomIds={invalidAtomIds}
          invalidBondIds={invalidBondIds}
        />

        <OrbitControls
          enablePan={true}
          enableZoom={true}
          enableRotate={true}
          minDistance={2}
          maxDistance={20}
        />
      </Canvas>
      <div className="absolute top-2 left-2 text-xs opacity-60 bg-neutral-900 px-2 py-1 rounded z-10">
        3D Viewer
      </div>
    </div>
  );
}
