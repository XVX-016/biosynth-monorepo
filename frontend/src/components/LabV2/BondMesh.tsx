import { useMemo } from "react";
import * as THREE from "three";

interface BondMeshProps {
    aPos: number[];
    bPos: number[];
    order?: number;
}

export default function BondMesh({ aPos, bPos, order = 1 }: BondMeshProps) {
    const start = useMemo(() => new THREE.Vector3(aPos[0], aPos[1], aPos[2]), [aPos]);
    const end = useMemo(() => new THREE.Vector3(bPos[0], bPos[1], bPos[2]), [bPos]);
    const diff = new THREE.Vector3().subVectors(end, start);
    const length = diff.length();

    // Midpoint
    const mid = new THREE.Vector3().addVectors(start, end).multiplyScalar(0.5);

    // Orientation quaternion from Y-axis to bond direction
    const orientation = useMemo(() => {
        const yAxis = new THREE.Vector3(0, 1, 0);
        const quaternion = new THREE.Quaternion();
        quaternion.setFromUnitVectors(yAxis, diff.clone().normalize());
        return quaternion;
    }, [diff]);

    // For double/triple bonds, calculate perpendicular offset
    const perpendicular = useMemo(() => {
        const dir = diff.clone().normalize();
        // Find a perpendicular vector
        const arbitrary = Math.abs(dir.y) < 0.9 ? new THREE.Vector3(0, 1, 0) : new THREE.Vector3(1, 0, 0);
        const perp = new THREE.Vector3().crossVectors(dir, arbitrary).normalize();
        return perp;
    }, [diff]);

    const bondRadius = 0.06;
    const bondOffset = 0.15; // Distance between parallel bonds

    if (order === 1) {
        // Single bond - one cylinder
        return (
            <mesh position={[mid.x, mid.y, mid.z]} quaternion={orientation}>
                <cylinderGeometry args={[bondRadius, bondRadius, length, 12]} />
                <meshStandardMaterial color="#888" metalness={0.2} roughness={0.5} />
            </mesh>
        );
    }

    if (order === 2) {
        // Double bond - two parallel cylinders
        const offset1 = perpendicular.clone().multiplyScalar(bondOffset);
        const offset2 = perpendicular.clone().multiplyScalar(-bondOffset);

        return (
            <>
                <mesh
                    position={[mid.x + offset1.x, mid.y + offset1.y, mid.z + offset1.z]}
                    quaternion={orientation}
                >
                    <cylinderGeometry args={[bondRadius * 0.8, bondRadius * 0.8, length, 12]} />
                    <meshStandardMaterial color="#888" metalness={0.2} roughness={0.5} />
                </mesh>
                <mesh
                    position={[mid.x + offset2.x, mid.y + offset2.y, mid.z + offset2.z]}
                    quaternion={orientation}
                >
                    <cylinderGeometry args={[bondRadius * 0.8, bondRadius * 0.8, length, 12]} />
                    <meshStandardMaterial color="#888" metalness={0.2} roughness={0.5} />
                </mesh>
            </>
        );
    }

    if (order === 3) {
        // Triple bond - three parallel cylinders
        const offset1 = perpendicular.clone().multiplyScalar(bondOffset);
        const offset2 = perpendicular.clone().multiplyScalar(-bondOffset);

        return (
            <>
                {/* Center cylinder */}
                <mesh position={[mid.x, mid.y, mid.z]} quaternion={orientation}>
                    <cylinderGeometry args={[bondRadius * 0.7, bondRadius * 0.7, length, 12]} />
                    <meshStandardMaterial color="#888" metalness={0.2} roughness={0.5} />
                </mesh>
                {/* Offset cylinders */}
                <mesh
                    position={[mid.x + offset1.x, mid.y + offset1.y, mid.z + offset1.z]}
                    quaternion={orientation}
                >
                    <cylinderGeometry args={[bondRadius * 0.7, bondRadius * 0.7, length, 12]} />
                    <meshStandardMaterial color="#888" metalness={0.2} roughness={0.5} />
                </mesh>
                <mesh
                    position={[mid.x + offset2.x, mid.y + offset2.y, mid.z + offset2.z]}
                    quaternion={orientation}
                >
                    <cylinderGeometry args={[bondRadius * 0.7, bondRadius * 0.7, length, 12]} />
                    <meshStandardMaterial color="#888" metalness={0.2} roughness={0.5} />
                </mesh>
            </>
        );
    }

    // Fallback for higher orders - just render as single thick bond
    return (
        <mesh position={[mid.x, mid.y, mid.z]} quaternion={orientation}>
            <cylinderGeometry args={[bondRadius * 1.2, bondRadius * 1.2, length, 12]} />
            <meshStandardMaterial color="#888" metalness={0.2} roughness={0.5} />
        </mesh>
    );
}
