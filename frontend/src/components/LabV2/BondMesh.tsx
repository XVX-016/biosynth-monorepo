import * as THREE from "three";

/** Render a cylinder between two atom positions */
export default function BondMesh({ aPos, bPos }: { aPos: number[]; bPos: number[] }) {
    const [ax, ay, az] = aPos;
    const [bx, by, bz] = bPos;

    // Compute midpoint and orientation
    const start = new THREE.Vector3(ax, ay, az);
    const end = new THREE.Vector3(bx, by, bz);
    const diff = new THREE.Vector3().subVectors(end, start);
    const length = diff.length();
    const midpoint = new THREE.Vector3().addVectors(start, end).multiplyScalar(0.5);

    // Create quaternion for rotation
    const direction = diff.clone().normalize();
    const axis = new THREE.Vector3(0, 1, 0);
    const quaternion = new THREE.Quaternion().setFromUnitVectors(axis, direction);

    return (
        <mesh position={midpoint} quaternion={quaternion}>
            <cylinderGeometry args={[0.06, 0.06, length, 10]} />
            <meshStandardMaterial color="#888" />
        </mesh>
    );
}
