import { Suspense, useRef } from 'react';
import { Canvas, useFrame } from '@react-three/fiber';
import { Environment, Float, PerspectiveCamera, OrbitControls } from '@react-three/drei';
import * as THREE from 'three';
import type { MoleculeGraph } from '@biosynth/engine';
import { moleculeToRenderable } from '../../lib/engineAdapter';

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

    useFrame((state) => {
        if (!groupRef.current) return;
        groupRef.current.rotation.y += 0.005;
        const time = state.clock.elapsedTime;
        groupRef.current.position.y = Math.sin(time * 0.5) * 0.1;
    });

    return (
        <group ref={groupRef} scale={[0.5, 0.5, 0.5]}>
            {renderable.atoms.map((atom) => {
                const radius = ELEMENT_RADII[atom.element] || 1.0;
                return (
                    <mesh key={atom.id} position={atom.position as any}>
                        <sphereGeometry args={[radius, 32, 32]} />
                        <meshPhysicalMaterial
                            color="#ffffff"
                            metalness={0.9}
                            roughness={0.1}
                            envMapIntensity={1.5}
                        />
                    </mesh>
                );
            })}

            {renderable.bonds.map((bond) => {
                const vFrom = new THREE.Vector3(...bond.from as [number, number, number]);
                const vTo = new THREE.Vector3(...bond.to as [number, number, number]);
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
                        <cylinderGeometry args={[radius, radius, length, 32]} />
                        <meshPhysicalMaterial
                            color="#ffffff"
                            metalness={0.9}
                            roughness={0.1}
                            envMapIntensity={1.5}
                        />
                    </mesh>
                );
            })}
        </group>
    );
}

interface MentorAvatarProps {
    modelPath: string;
    accentColor: string;
}

function MentorAvatar({ accentColor }: MentorAvatarProps) {
    const groupRef = useRef<THREE.Group>(null);
    const headRef = useRef<THREE.Mesh>(null);

    // Simple head-tracking logic
    useFrame((state) => {
        if (!headRef.current) return;

        // Target position based on mouse
        const targetX = (state.mouse.x * Math.PI) / 8;
        const targetY = (state.mouse.y * Math.PI) / 10;

        // Smooth Interpolation
        headRef.current.rotation.y = THREE.MathUtils.lerp(headRef.current.rotation.y, targetX, 0.1);
        headRef.current.rotation.x = THREE.MathUtils.lerp(headRef.current.rotation.x, -targetY, 0.1);

        if (groupRef.current) {
            // Idle breathing
            groupRef.current.position.y = Math.sin(state.clock.elapsedTime * 0.8) * 0.05 - 1.2;
        }
    });

    return (
        <group ref={groupRef} position={[0, -1.2, 0]}>
            {/* Placeholder for Upper-Body Avatar */}
            {/* Chest/Shoulders */}
            <mesh position={[0, 0, 0]}>
                <capsuleGeometry args={[0.5, 0.8, 4, 16]} />
                <meshPhysicalMaterial
                    color={accentColor}
                    metalness={0.1}
                    roughness={0.8}
                    transmission={0.5}
                    thickness={1}
                    opacity={0.3}
                    transparent
                />
            </mesh>

            {/* Head */}
            <mesh ref={headRef} position={[0, 0.8, 0]}>
                <sphereGeometry args={[0.3, 32, 32]} />
                <meshPhysicalMaterial
                    color="#ffffff"
                    metalness={0.2}
                    roughness={0.1}
                    envMapIntensity={2}
                />
                {/* "Eyes" for directionality */}
                <mesh position={[0.1, 0.05, 0.25]}>
                    <sphereGeometry args={[0.04, 16, 16]} />
                    <meshBasicMaterial color="#000000" />
                </mesh>
                <mesh position={[-0.1, 0.05, 0.25]}>
                    <sphereGeometry args={[0.04, 16, 16]} />
                    <meshBasicMaterial color="#000000" />
                </mesh>
            </mesh>
        </group>
    );
}

interface Studio3DSceneProps {
    mentorId: string;
    accentColor: string;
    modelPath: string;
    molecule?: MoleculeGraph | null;
}

export default function Studio3DScene({ accentColor, modelPath, molecule }: Studio3DSceneProps) {
    return (
        <div className="w-full h-full relative">
            <Canvas shadows dpr={[1, 2]}>
                <PerspectiveCamera makeDefault position={[0, 0, 4]} fov={45} />

                <Suspense fallback={null}>
                    <Environment preset="city" />
                    <ambientLight intensity={0.5} />
                    <pointLight position={[10, 10, 10]} intensity={1} castShadow />

                    <Float speed={1.5} rotationIntensity={0.5} floatIntensity={0.5}>
                        <MentorAvatar accentColor={accentColor} modelPath={modelPath} />
                    </Float>

                    {molecule && (
                        <group position={[1.4, 0.4, 0.5]}>
                            <FloatingMolecule molecule={molecule} />
                        </group>
                    )}

                    <OrbitControls
                        enableZoom={false}
                        enablePan={false}
                        autoRotate={!!molecule}
                        autoRotateSpeed={0.5}
                    />

                    {/* Background Glow */}
                    <mesh position={[0, 0, -2]}>
                        <planeGeometry args={[10, 10]} />
                        <meshBasicMaterial
                            color={accentColor}
                            transparent
                            opacity={0.05}
                            side={THREE.DoubleSide}
                        />
                    </mesh>
                </Suspense>
            </Canvas>
        </div>
    );
}
