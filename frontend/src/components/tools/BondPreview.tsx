import React, { useRef } from 'react';
import { useFrame, useThree } from '@react-three/fiber';
import * as THREE from 'three';
import { useMoleculeStore } from '../../store/moleculeStore';

export function BondPreview() {
    const { camera, mouse, scene } = useThree();
    const lineRef = useRef<any>(null);
    const selectedAtomId = useMoleculeStore((state) => state.selectedAtomId);
    const tool = useMoleculeStore((state) => state.tool);
    const currentMolecule = useMoleculeStore((state) => state.currentMolecule);

    useFrame(() => {
        if (!lineRef.current || tool !== "bond" || !selectedAtomId || !currentMolecule) {
            if (lineRef.current) lineRef.current.visible = false;
            return;
        }

        const atom = currentMolecule.atoms.get(selectedAtomId);
        if (!atom) return;

        lineRef.current.visible = true;

        // Start point: Atom position
        const start = new THREE.Vector3(...atom.position);

        // End point: Mouse raycast to plane (similar to add-atom logic)
        // Note: Ideally we raycast to other atoms, but for loose space, a plane at z=atom.z is fine or camera plane
        // Let's simple unproject
        const vector = new THREE.Vector3(mouse.x, mouse.y, 0.5);
        vector.unproject(camera);
        const dir = vector.sub(camera.position).normalize();
        const distance = -camera.position.z / dir.z; // Project to Z=0 plane default?
        // Better: Project to plane parallel to view at atom depth? 
        // Or just some distance in front of camera

        // Simple approach: Raycast to Z=0 plane (matches LabCanvas logic)
        // If atom is not at Z=0, this looks weird. But better than nothing.
        const planeZ = atom.position[2]; // Project to same depth as start atom
        // Ray-plane intersection
        const t = (planeZ - camera.position.z) / dir.z;
        const pos = camera.position.clone().add(dir.multiplyScalar(t));

        const positions = new Float32Array([
            start.x, start.y, start.z,
            pos.x, pos.y, pos.z
        ]);

        lineRef.current.geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    });

    if (tool !== "bond" || !selectedAtomId) return null;

    return (
        <line ref={lineRef}>
            <bufferGeometry />
            <lineBasicMaterial color="#00ff00" transparent opacity={0.5} />
        </line>
    );
}
