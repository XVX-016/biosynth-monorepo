import React from "react";
import { apiClient } from "../../lib/api";
import { useEditor } from "../../context/EditorContext";
import { useNavigate } from "react-router-dom";

export default function ImportToLabButton({ moleculeId }: { moleculeId: string }) {
    const { dispatch } = useEditor();
    const navigate = useNavigate();

    const handleImport = async () => {
        try {
            // Use apiClient which has base URL configured
            const res = await apiClient.post(`/library/molecules/${moleculeId}/export`);
            const mol = res.data.data ?? res.data;

            if (mol.atoms) {
                const atoms = mol.atoms.map((a: any, i: number) => ({
                    id: a.id || `atom_${i}`,
                    element: a.element,
                    position: [a.x, a.y, a.z || 0]
                }));
                const bonds = (mol.bonds || []).map((b: any, i: number) => ({
                    id: b.id || `bond_${i}`,
                    a: b.a,
                    b: b.b,
                    order: b.order
                }));
                dispatch({ type: "LOAD_MOLECULE", payload: { atoms, bonds, name: mol.name } });
                navigate("/lab");
            }
        } catch (e) {
            console.error("Import failed", e);
        }
    };

    return <button className="px-2 py-1 bg-gray-900 text-white rounded text-sm hover:bg-black" onClick={handleImport}>Open in Lab</button>;
}
