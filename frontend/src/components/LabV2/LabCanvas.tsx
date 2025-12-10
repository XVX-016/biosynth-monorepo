import { useCallback, Suspense } from "react";
import { Canvas } from "@react-three/fiber";
import { OrbitControls, Html, Grid } from "@react-three/drei";
import { useLabStore } from "../../store/labStore";
import AtomMesh from "./AtomMesh";
import BondMesh from "./BondMesh";

function Scene() {
    const { molecule, currentTool, addAtom } = useLabStore();

    // Click on plane to add atom (only when tool = add-atom)
    const onPlanePointerDown = useCallback((e: any) => {
        if (currentTool !== "add-atom" && currentTool !== "element") return;
        e.stopPropagation();
        const point = e.point;
        addAtom("C", [point.x, point.y, point.z]);
    }, [currentTool, addAtom]);

    // Use Drei Grid with light styling to match user request
    return (
        <>
            <ambientLight intensity={0.9} />
            <directionalLight position={[5, 10, 7]} intensity={0.6} />

            {/* Invisible plane for pointer placement */}
            <mesh rotation={[-Math.PI / 2, 0, 0]} position={[0, -0.01, 0]} receiveShadow onPointerDown={onPlanePointerDown}>
                <planeGeometry args={[100, 100]} />
                <meshBasicMaterial transparent opacity={0} />
            </mesh>

            {/* Light Grid */}
            <Grid args={[20, 20]} cellColor="#e5e7eb" sectionColor="#cbd5e1" fadeDistance={30} infiniteGrid />

            {/* Bonds */}
            {molecule.bonds.map((b) => {
                const a = molecule.atoms.find((x) => x.id === b.atom1);
                const c = molecule.atoms.find((x) => x.id === b.atom2);
                if (!a || !c) return null;
                return <BondMesh key={b.id} aPos={a.position} bPos={c.position} order={b.order} />;
            })}

            {/* Atoms */}
            {molecule.atoms.map(a => <AtomMesh key={a.id} atom={a} />)}
        </>
    );
}

export default function LabCanvas() {
    return (
        <Canvas camera={{ position: [0, 6, 10], fov: 50 }} style={{ width: "100%", height: "100%" }} gl={{ alpha: true }}>
            <OrbitControls makeDefault />
            <Suspense fallback={<Html center><div className="text-gray-500">Loading...</div></Html>}>
                <Scene />
            </Suspense>
        </Canvas>
    );
}
