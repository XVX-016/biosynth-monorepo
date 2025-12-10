import React from "react";
import { useEditor } from "../../../context/EditorContext";

export default function AddAtomButton() {
    const { dispatch, state } = useEditor();
    const isActive = state.tool === "add-atom";
    return (
        <button
            className={`px-2 py-1 border rounded ${isActive ? "bg-blue-100 border-blue-500" : "hover:bg-gray-50"}`}
            onClick={() => dispatch({ type: "SET_TOOL", payload: "add-atom" })}
        >
            Add Atom
        </button>
    );
}
