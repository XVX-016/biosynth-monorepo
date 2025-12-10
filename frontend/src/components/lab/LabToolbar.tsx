import React from "react";
import OptimizeButton from "./OptimizeButton";
import PredictButton from "./PredictButton";
import AutoBondButton from "./tools/AutoBondButton";
import AddAtomButton from "./tools/AddAtomButton";
import AddBondButton from "./tools/AddBondButton";
import ExportButton from "./ExportButton";
import { useEditor } from "../../context/EditorContext";

// Stub OptimizeButton since we didn't create it explicitly yet
const OptimizeButtonStub = () => <button className="px-2 py-1 border rounded hover:bg-gray-50 text-gray-400 cursor-not-allowed">Optimize (Stub)</button>;

export default function LabToolbar() {
    const { undo, redo, canUndo, canRedo } = useEditor();

    return (
        <div className="w-full border-b border-gray-300 bg-white p-2 flex items-center gap-2 overflow-x-auto">
            <AddAtomButton />
            <AddBondButton />
            <AutoBondButton />

            <div className="w-px h-6 bg-gray-300 mx-2" />

            <button onClick={() => undo()} disabled={!canUndo()} className="px-2 py-1 border rounded disabled:opacity-50 hover:bg-gray-50">Undo</button>
            <button onClick={() => redo()} disabled={!canRedo()} className="px-2 py-1 border rounded disabled:opacity-50 hover:bg-gray-50">Redo</button>

            <div className="w-px h-6 bg-gray-300 mx-2" />

            <OptimizeButton />
            <PredictButton />
            <ExportButton />
        </div>
    );
}
