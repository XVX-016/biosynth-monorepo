import React from "react";
import { useEditorContext } from "../../../context/EditorContext";

export default function AutoBondButton() {
    const { state, dispatch } = useEditorContext();
    return (
        <button
            className={`btn ${state.autoBond ? "btn-active" : ""}`}
            onClick={() => dispatch({ type: "TOGGLE_AUTOBOND" })}
        >
            AutoBond
        </button>
    );
}
