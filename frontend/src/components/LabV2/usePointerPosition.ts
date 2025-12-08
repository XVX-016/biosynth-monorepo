import { useThree } from "@react-three/fiber";
import * as THREE from "three";

/**
 * Returns getPointerOnPlane(nativeMouseEvent) -> {x,y,z}
 * Projects pointer onto y=0 plane (approximate interaction plane)
 */
export default function usePointerPosition() {
    const { camera, gl } = useThree();

    function getPointerOnPlane(e: MouseEvent) {
        const rect = (gl.domElement as HTMLCanvasElement).getBoundingClientRect();
        const x = ((e.clientX - rect.left) / rect.width) * 2 - 1;
        const y = -((e.clientY - rect.top) / rect.height) * 2 + 1;

        // Raycast
        const v = new THREE.Vector3(x, y, 0.5).unproject(camera);
        const dir = v.sub(camera.position).normalize();
        const distance = -camera.position.y / dir.y; // intersect y=0 plane

        if (!isFinite(distance) || distance < 0) return null;

        const pos = camera.position.clone().add(dir.multiplyScalar(distance));
        return { x: pos.x, y: 0, z: pos.z };
    }

    return { getPointerOnPlane };
}
