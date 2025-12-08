import { Suspense, useRef, useEffect } from "react";
import { Canvas } from "@react-three/fiber";
import { OrbitControls, Html, Grid } from "@react-three/drei";
import AtomMesh from "./AtomMesh";
import BondMesh from "./BondMesh";
import { useEditor } from "../../context/EditorContext";
import { nanoid } from "nanoid";

/**
 * LabCanvas is the 3D view.
 * - renders atoms and bonds from EditorContext
 * - allows clicking to add atoms when tool='add-atom'
 * - shows grid and subtle lighting
 */

function SceneContents() {
    const { state, dispatch } = useEditor();
    const controlsRef = useRef<any>(null);
    const planeRef = useRef<any>(null);

    // Auto-fit camera when molecule changes
    useEffect(() => {
        if (state.atoms && state.atoms.length > 0) {
            // Simple camera positioning based on atom bounds
            const positions = state.atoms.map((a: any) => {
                const pos = a.position || [a.x || 0, a.y || 0, a.z || 0];
                return pos;
            });

            // Calculate center
            const center = positions.reduce(
                (acc, pos) => [acc[0] + pos[0], acc[1] + pos[1], acc[2] + pos[2]],
                [0, 0, 0]
            ).map((v: number) => v / positions.length);

            if (controlsRef.current) {
                controlsRef.current.target.set(center[0], center[1], center[2]);
            }
        }
    }, [state.atoms?.length]);

    const handlePlaneClick = (e: any) => {
        if (state.tool !== "add-atom") return;

        e.stopPropagation();

        // Get world position from intersection
        const point = e.point;

        // Create atom at world coords
        const atom = {
            id: nanoid(),
            element: "C",
            position: [point.x, point.y, point.z] as [number, number, number]
        };

        dispatch({ type: "ADD_ATOM", payload: atom });

        // If auto bond on, trigger autobond
        if (state.autoBond) {
            dispatch({ type: "AUTO_BOND" });
        }
    };

    return (
        <>
            {/* Lighting */}
            <ambientLight intensity={0.6} />
            <directionalLight position={[5, 10, 7]} intensity={0.6} />

            {/* Grid floor */}
            <Grid
                args={[40, 40]}
                cellColor="#333"
                sectionColor="#444"
                fadeDistance={50}
                fadeStrength={1}
                position={[0, -0.01, 0]}
            />

            {/* Invisible interaction plane */}
            <mesh
                ref={planeRef}
                rotation={[-Math.PI / 2, 0, 0]}
                position={[0, 0, 0]}
                onClick={handlePlaneClick}
            >
                <planeGeometry args={[100, 100]} />
                <meshBasicMaterial transparent opacity={0} />
            </mesh>

            {/* Render bonds */}
            {state.bonds?.map((b: any) => {
                const a = state.atoms?.find((x: any) => x.id === b.a);
                const c = state.atoms?.find((x: any) => x.id === b.b);
                if (!a || !c) return null;

                const aPos = a.position || [a.x || 0, a.y || 0, a.z || 0];
                const bPos = c.position || [c.x || 0, c.y || 0, c.z || 0];

                return <BondMesh key={b.id ?? `${b.a}_${b.b}`} aPos={aPos} bPos={bPos} />;
            })}

            {/* Render atoms */}
            {state.atoms?.map((a: any) => (
                <AtomMesh key={a.id} atom={a} />
            ))}

            <OrbitControls ref={controlsRef} makeDefault />
        </>
    );
}

export default function LabCanvas() {
    return (
        <Canvas camera={{ position: [0, 6, 10], fov: 50 }} style={{ background: "#0a0a0a" }}>
            <Suspense fallback={<Html center><div style={{ color: "white" }}>Loading...</div></Html>}>
                <SceneContents />
            </Suspense>
        </Canvas>
    );
}
