import React, { Suspense } from 'react';
import { Canvas } from '@react-three/fiber';
import { OrbitControls } from '@react-three/drei';
import AtomMesh from './r3f/AtomMesh';
import BondMesh from './r3f/BondMesh';

interface MoleculeCanvasProps {
    atoms: any[];
    bonds: any[];
}

export default function MoleculeCanvas({ atoms, bonds }: MoleculeCanvasProps) {
    if (!atoms || !bonds) return null;

    return (
        <div className="w-full h-full bg-gray-50">
            <Canvas camera={{ position: [0, 0, 15], fov: 45 }}>
                <ambientLight intensity={0.6} />
                <directionalLight position={[10, 10, 10]} intensity={1} />
                <Suspense fallback={null}>
                    {atoms.map((atom: any) => (
                        <AtomMesh
                            key={atom.id}
                            id={atom.id}
                            position={atom.position}
                            element={atom.element}
                            renderMode="ballstick"
                            colorScheme="element"
                        />
                    ))}
                    {bonds.map((bond: any) => (
                        <BondMesh
                            key={bond.id}
                            id={bond.id}
                            from={bond.a1 || bond.a} // Handle both formats
                            to={bond.a2 || bond.b}
                            order={bond.order}
                            renderMode="ballstick"
                        />
                    ))}
                    <OrbitControls makeDefault />
                </Suspense>
            </Canvas>
        </div>
    );
}
