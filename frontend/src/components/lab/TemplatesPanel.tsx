import React from "react";
import { CH4, BENZENE } from "../../utils/defaults";
import { useEditorContext } from "../../context/EditorContext";

export default function TemplatesPanel() {
    const { dispatch } = useEditorContext();

    const insert = (mol: any) => {
        // ensure unique ids by prefixing
        const prefix = Date.now().toString().slice(-4);
        const atoms = mol.atoms.map((a: any) => ({ ...a, id: `${prefix}_${a.id}` }));
        const bonds = mol.bonds.map((b: any) => ({ ...b, a: `${prefix}_${b.a}`, b: `${prefix}_${b.b}` }));
        dispatch({ type: "SET_MOLECULE", payload: { atoms, bonds, name: mol.name } });
    };

    return (
        <div className="p-2">
            <div className="font-bold mb-2">Templates</div>
            <button className="btn block mb-1" onClick={() => insert(CH4)}>Methane (CH4)</button>
            <button className="btn block" onClick={() => insert(BENZENE)}>Benzene</button>
        </div>
    );
}
