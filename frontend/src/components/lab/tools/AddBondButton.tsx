import React from "react";
import { useEditorContext } from "../../../context/EditorContext";

export default function AddBondButton() {
    const { dispatch } = useEditorContext();
    return (
        <button className="btn" onClick={() => dispatch({ type: "SET_TOOL", payload: "add-bond" })}>
            Add Bond
        </button>
    );
}
