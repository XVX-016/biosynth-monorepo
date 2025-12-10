import React, { useEffect, useRef } from 'react';
import * as THREE from 'three';
import Renderer3D from '../../lab/renderer/Renderer3D';
import { useEditor } from '../../context/EditorContext';

export default function Lab3DViewport() {
    const mount = useRef<HTMLDivElement | null>(null);
    const rendererRef = useRef<Renderer3D | null>(null);
    const { state } = useEditor();

    useEffect(() => {
        if (!mount.current) return;
        const r = new Renderer3D(mount.current);
        rendererRef.current = r;
        r.start();

        // attach basic pointer handler
        const onClick = (e: PointerEvent) => r.onPointerDown(e);
        mount.current.addEventListener('pointerdown', onClick as any);

        return () => {
            mount.current?.removeEventListener('pointerdown', onClick as any);
            r.dispose();
        };
    }, []);

    // Sync atoms: add missing atom spheres
    useEffect(() => {
        const r = rendererRef.current;
        if (!r) return;
        // add atoms not yet present
        const existing = new Set(Array.from(r.atomsMap.keys()));
        for (const a of state.atoms) {
            if (!existing.has(a.id)) {
                r.addAtomSphere(a.id, a.element, a.position);
            } else {
                // update position if needed
                const obj = r.atomsMap.get(a.id) as THREE.Object3D;
                if (obj) obj.position.set(...a.position);
            }
        }
        // remove atoms that no longer exist
        for (const id of Array.from(r.atomsMap.keys())) {
            if (!state.atoms.find(x => x.id === id)) {
                const obj = r.atomsMap.get(id);
                if (obj) r.scene.remove(obj);
                r.atomsMap.delete(id);
            }
        }
    }, [state.atoms]);

    // Sync bonds: animate new bonds
    useEffect(() => {
        const r = rendererRef.current;
        if (!r) return;
        const existing = new Set(Array.from(r.bondsMap.keys()));
        for (const b of state.bonds) {
            if (!existing.has(b.id)) {
                // find atom positions
                const aObj = r.atomsMap.get(b.a) as THREE.Object3D | undefined;
                const bObj = r.atomsMap.get(b.b) as THREE.Object3D | undefined;
                if (aObj && bObj) {
                    const posA: [number, number, number] = [aObj.position.x, aObj.position.y, aObj.position.z];
                    const posB: [number, number, number] = [bObj.position.x, bObj.position.y, bObj.position.z];
                    // animated creation
                    r.addBondAnimated(b.id, posA, posB, b.order).catch(console.error);
                }
            }
        }
        // remove bonds that no longer exist in state
        for (const id of Array.from(r.bondsMap.keys())) {
            if (!state.bonds.find(x => x.id === id)) {
                r.removeBond(id);
            }
        }
    }, [state.bonds]);

    return <div ref={mount} style={{ width: '100%', height: '100%' }} />;
}
