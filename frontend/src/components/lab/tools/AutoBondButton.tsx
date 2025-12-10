import React from "react";
import { useEditor } from "../../../context/EditorContext";

export default function AutoBondButton() {
    const { state, dispatch } = useEditor();
    return (
        <button
            className={`px-2 py-1 border rounded ${state.autoBond ? "bg-green-100 border-green-500 text-green-700" : "hover:bg-gray-50"}`}
            onClick={() => dispatch({ type: "TOGGLE_AUTOBOND" })}
        >
            AutoBond
        </button>
    );
}
