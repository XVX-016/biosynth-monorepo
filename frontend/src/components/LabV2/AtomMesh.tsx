import { useRef, useState, useCallback } from "react";
import { useLabStore } from "../../store/labStore";
import { getAtomColor } from "../../utils/atomColors";
import { getElementSpec } from "../../utils/elements";

export default function AtomMesh({ atom }: { atom: any }) {
    const mesh = useRef<any>();
    const { currentTool, addBond, deleteAtom, setSelectedAtomId, selectedAtomId, bondOrder } = useLabStore();
    const [hovered, setHovered] = useState(false);

    const color = getAtomColor(atom.element);
    const spec = getElementSpec(atom.element);
    // Use element radius if available, otherwise fallback to size based on element
    const radius = spec ? spec.radius * 0.25 : (atom.element === "H" ? 0.18 : 0.28);

    const handleClick = useCallback((e: any) => {
        e.stopPropagation();

        if (currentTool === "add-bond") {
            // Check if another atom is selected
            if (selectedAtomId && selectedAtomId !== atom.id) {
                // Create bond with selected bond order
                addBond(selectedAtomId, atom.id, bondOrder);
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
    }, [currentTool, selectedAtomId, atom.id, addBond, deleteAtom, setSelectedAtomId, bondOrder]);

    // Position from object (Core Atom)
    const { x, y, z } = atom.position;

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
