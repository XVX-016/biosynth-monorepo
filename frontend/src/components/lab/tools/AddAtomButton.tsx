import React from "react";
import { useEditorContext } from "../../../context/EditorContext";

export default function AddAtomButton() {
    const { dispatch } = useEditorContext();
    return (
        <button className="btn" onClick={() => dispatch({ type: "SET_TOOL", payload: "add-atom" })}>
            Add Atom
        </button>
    );
}
