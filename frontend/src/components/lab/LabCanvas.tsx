import { Canvas } from "@react-three/fiber";
import { OrbitControls } from "@react-three/drei";

export function LabCanvas() {
    return (
        <Canvas
            camera={{ position: [5, 5, 5], fov: 45 }}
            shadows
            style={{ width: "100%", height: "100%" }}
        >
            {/* Safe Lighting */}
            <ambientLight intensity={0.8} />
            <directionalLight intensity={1} position={[5, 10, 5]} castShadow />

            {/* White Plane */}
            <mesh rotation={[-Math.PI / 2, 0, 0]} receiveShadow>
                <planeGeometry args={[100, 100]} />
                <meshStandardMaterial color="#ffffff" />
            </mesh>

            <OrbitControls enableDamping />
        </Canvas>
    );
}
