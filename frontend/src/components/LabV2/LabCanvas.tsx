import { useCallback, Suspense, useMemo } from "react";
import { Canvas } from "@react-three/fiber";
import { OrbitControls, Html, Line } from "@react-three/drei";
import { useLabStore } from "../../store/labStore";
import AtomMesh from "./AtomMesh";
import BondMesh from "./BondMesh";

function CustomGrid() {
    const lines = useMemo(() => {
        const items = [];
        const size = 100;     // Total width
        const spacing = 2;    // Spacing between lines
        const steps = size / spacing;

        // Draw lines parallel to X-axis (Horizontal)
        // These vary in Z
        for (let i = -steps / 2; i <= steps / 2; i++) {
            const z = i * spacing;
            items.push(
                <Line
                    key={`h-${i}`}
                    points={[[-size / 2, 0, z], [size / 2, 0, z]]}
                    color="#e5e7eb" // Tailwind gray-200
                    lineWidth={1}
                    transparent
                    opacity={0.4}
                />
            );
        }

        // Draw lines parallel to Z-axis (Vertical)
        // These vary in X
        for (let i = -steps / 2; i <= steps / 2; i++) {
            const x = i * spacing;
            items.push(
                <Line
                    key={`v-${i}`}
                    points={[[x, 0, -size / 2], [x, 0, size / 2]]}
                    color="#e5e7eb" // Tailwind gray-200
                    lineWidth={1}
                    transparent
                    opacity={0.4}
                />
            );
        }
        return items;
    }, []);

    return <group position={[0, -0.01, 0]}>{lines}</group>;
}

function Scene() {
    const { molecule, currentTool, currentElement, addAtom } = useLabStore();

    // Click on plane to add atom (only when tool = add-atom)
    const onPlanePointerDown = useCallback((e: any) => {
        if (currentTool !== "add-atom") return;
        e.stopPropagation();
        const point = e.point;
        // Basic snapping could be added here
        addAtom(currentElement, [point.x, point.y, point.z]);
    }, [currentTool, addAtom, currentElement]);

    return (
        <>
            <ambientLight intensity={0.9} />
            <directionalLight position={[5, 10, 7]} intensity={0.6} />
            <directionalLight position={[-5, 5, -5]} intensity={0.3} />

            {/* Invisible plane for pointer interaction */}
            <mesh rotation={[-Math.PI / 2, 0, 0]} position={[0, -0.01, 0]} receiveShadow onPointerDown={onPlanePointerDown}>
                <planeGeometry args={[100, 100]} />
                <meshBasicMaterial transparent opacity={0} />
            </mesh>

            {/* Custom Grid with both Horizontal and Vertical lines */}
            <CustomGrid />

            {/* Center Origin Dot */}
            <mesh position={[0, 0, 0]} rotation={[-Math.PI / 2, 0, 0]}>
                <circleGeometry args={[0.05, 32]} />
                <meshBasicMaterial color="#9ca3af" opacity={0.5} transparent />
            </mesh>

            {/* Bonds */}
            {molecule.bonds.map((b) => {
                const a = molecule.atoms.find((x) => x.id === b.from);
                const c = molecule.atoms.find((x) => x.id === b.to);
                if (!a || !c) return null;
                // Atoms have position: {x,y,z}. Map to array for BondMesh if needed?
                // BondMesh expects [x,y,z] arrays?
                // Let's check BondMesh props.
                // Assuming BondMesh updated or we map.
                // Core Atom has .position as object. App Molecule extends Core.
                // BondMesh likely expects [n,n,n].
                // We need to pass [a.position.x, a.position.y, a.position.z].

                return <BondMesh key={b.id}
                    aPos={[a.position!.x, a.position!.y, a.position!.z]}
                    bPos={[c.position!.x, c.position!.y, c.position!.z]}
                    order={b.order}
                />;
            })}

            {/* Atoms */}
            {molecule.atoms.map(a => <AtomMesh key={a.id} atom={a} />)}
        </>
    );
}

export default function LabCanvas() {
    return (
        <Canvas
            camera={{ position: [0, 6, 12], fov: 45 }}
            style={{ width: "100%", height: "100%" }}
            gl={{ alpha: true, antialias: true }}
            dpr={[1, 2]} // Handle high DPI
        >
            <OrbitControls
                makeDefault
                enableDamping
                dampingFactor={0.1}
                maxPolarAngle={Math.PI / 2 - 0.1} // Don't go below ground
                minDistance={2}
                maxDistance={50}
            />
            <Scene />
        </Canvas>
    );
}
