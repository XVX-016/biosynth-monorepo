import React from "react";
import { CH4, BENZENE } from "../../utils/defaultMolecules";
import { useEditor } from "../../context/EditorContext";

export default function TemplatesPanel() {
    const { dispatch } = useEditor();

    const insert = (mol: any) => {
        const prefix = Date.now().toString().slice(-4);
        // map to state Atom format: position array
        const atoms = mol.atoms.map((a: any) => ({
            id: `${prefix}_${a.id}`,
            element: a.element,
            position: [a.x, a.y, a.z]
        }));
        const bonds = mol.bonds.map((b: any) => ({ ...b, id: `${prefix}_b_${b.a}_${b.b}`, a: `${prefix}_${b.a}`, b: `${prefix}_${b.b}` }));
        dispatch({ type: "LOAD_MOLECULE", payload: { atoms, bonds } });
    };

    return (
        <div className="flex flex-col gap-2">
            <h3 className="font-semibold text-sm">Templates</h3>
            <button className="px-2 py-1 text-sm border rounded hover:bg-gray-50 text-left" onClick={() => insert(CH4)}>Methane (CH4)</button>
            <button className="px-2 py-1 text-sm border rounded hover:bg-gray-50 text-left" onClick={() => insert(BENZENE)}>Benzene</button>
        </div>
    );
}
