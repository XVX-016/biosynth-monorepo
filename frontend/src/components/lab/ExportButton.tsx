import React from "react";
import { useEditorContext } from "../../context/EditorContext";
import { exportAsMol } from "../../utils/serializeMolecule";

export default function ExportButton() {
    const { state } = useEditorContext();
    const handleExport = () => {
        const mol = { atoms: state.atoms, bonds: state.bonds, name: state.name || "export" };
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

    return <button className="btn" onClick={handleExport}>Export .mol</button>;
}
