import React from "react";
import { useEditor } from "../../../context/EditorContext";

export default function AddBondButton() {
    const { dispatch, state } = useEditor();
    const isActive = state.tool === "add-bond";
    return (
        <button
            className={`px-2 py-1 border rounded ${isActive ? "bg-blue-100 border-blue-500" : "hover:bg-gray-50"}`}
            onClick={() => dispatch({ type: "SET_TOOL", payload: "add-bond" })}
        >
            Add Bond
        </button>
    );
}
