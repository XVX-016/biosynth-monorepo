import React, { useMemo } from 'react';
import * as THREE from 'three';

interface BondProps {
    start: [number, number, number];
    end: [number, number, number];
    order?: number;
    color?: string;
}

export function Bond({ start, end, order = 1, color = "gray" }: BondProps) {
    const { position, rotation, height } = useMemo(() => {
        const startVec = new THREE.Vector3(...start);
        const endVec = new THREE.Vector3(...end);

        const distance = startVec.distanceTo(endVec);
        const position = startVec.clone().add(endVec).multiplyScalar(0.5);

        const quaternion = new THREE.Quaternion();
        const cylinderUp = new THREE.Vector3(0, 1, 0);
        const direction = endVec.clone().sub(startVec).normalize();

        quaternion.setFromUnitVectors(cylinderUp, direction);
        const rotation = new THREE.Euler().setFromQuaternion(quaternion);

        return { position, rotation, height: distance };
    }, [start, end]);

    const radius = 0.05 + (order - 1) * 0.03; // Thicker for higher order

    return (
        <mesh position={position} rotation={rotation}>
            <cylinderGeometry args={[radius, radius, height, 8]} />
            <meshStandardMaterial color={color} />
        </mesh>
    );
}

export default Bond;
