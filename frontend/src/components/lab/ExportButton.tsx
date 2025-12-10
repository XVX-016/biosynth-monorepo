import React from "react";
import { useEditor } from "../../context/EditorContext";
import { exportAsMol } from "../../utils/serializeMolecule";

export default function ExportButton() {
    const { state } = useEditor();
    const handleExport = () => {
        // map state atom positions to x,y,z expected by exportAsMol
        const molAtoms = state.atoms.map(a => ({
            id: a.id,
            element: a.element,
            x: a.position[0],
            y: a.position[1],
            z: a.position[2]
        }));
        const mol = { atoms: molAtoms, bonds: state.bonds, name: state.name || "export" };
        const s = exportAsMol(mol);
        const blob = new Blob([s], { type: "chemical/x-mdl-molfile" });
        const url = URL.createObjectURL(blob);
        const a = document.createElement("a");
        a.href = url;
        a.download = `${mol.name || "molecule"}.mol`;
        document.body.appendChild(a);
        a.click();
        a.remove();
        URL.revokeObjectURL(url);
    };

    return <button className="px-2 py-1 border rounded hover:bg-gray-50" onClick={handleExport}>Export .mol</button>;
}
