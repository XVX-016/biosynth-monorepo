import { useRef, useState, useCallback } from "react";
import { useLabStore } from "../../store/labStore";

export default function AtomMesh({ atom }: { atom: any }) {
    const mesh = useRef<any>();
    const { currentTool, addBond, deleteAtom, setSelectedAtomId, selectedAtomId } = useLabStore();
    const [hovered, setHovered] = useState(false);

    const colorMap: Record<string, string> = {
        C: "#222", H: "#eee", O: "#e44", N: "#3aa", S: "#fc0", P: "#f90", Cl: "#2ca", F: "#9c6", Br: "#a52"
    };
    const color = colorMap[atom.element] ?? "#888";
    const radius = atom.element === "H" ? 0.18 : 0.28;

    const handleClick = useCallback((e: any) => {
        e.stopPropagation();

        if (currentTool === "add-bond") {
            // Check if another atom is selected
            if (selectedAtomId && selectedAtomId !== atom.id) {
                // Create bond
                addBond(selectedAtomId, atom.id, 1);
                setSelectedAtomId(null);
            } else {
                // Select this atom to start bonding
                setSelectedAtomId(atom.id);
            }
        } else if (currentTool === "erase") {
            deleteAtom(atom.id);
        } else {
            // Select
            setSelectedAtomId(atom.id);
        }
    }, [currentTool, selectedAtomId, atom.id, addBond, deleteAtom, setSelectedAtomId]);

    // Position from array
    const [x, y, z] = atom.position;

    return (
        <mesh
            ref={mesh}
            position={[x, y, z]}
            onPointerOver={(e) => { e.stopPropagation(); setHovered(true); }}
            onPointerOut={(e) => { e.stopPropagation(); setHovered(false); }}
            onClick={handleClick}
            scale={hovered || selectedAtomId === atom.id ? 1.2 : 1}
        >
            <sphereGeometry args={[radius, 32, 32]} />
            <meshStandardMaterial color={color} metalness={0.1} roughness={0.4} />
            {/* Hover/Selection halo */}
            {(hovered || selectedAtomId === atom.id) && (
                <mesh>
                    <sphereGeometry args={[radius * 1.22, 16, 16]} />
                    <meshBasicMaterial color={"#3b82f6"} opacity={0.3} transparent />
                </mesh>
            )}
        </mesh>
    );
}
