import * as THREE from "three";

/** Render a sphere for an atom with element-based coloring */
export default function AtomMesh({ atom }: { atom: any }) {
    const colorMap: Record<string, string> = {
        C: "#222",
        H: "#eee",
        O: "#e64",
        N: "#48a",
        S: "#fc0",
        P: "#f90",
        Cl: "#2ca"
    };
    const color = colorMap[atom.element] ?? "#999";
    const radius = atom.element === "H" ? 0.18 : 0.28;

    const position = atom.position || [atom.x || 0, atom.y || 0, atom.z || 0];

    return (
        <mesh position={new THREE.Vector3(...position)}>
            <sphereGeometry args={[radius, 24, 24]} />
            <meshStandardMaterial metalness={0.1} roughness={0.5} color={color} />
        </mesh>
    );
}
