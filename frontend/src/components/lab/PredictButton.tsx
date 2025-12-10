import React from "react";
import { useEditor } from "../../context/EditorContext";
import { MLAPI } from "../../api/ml";

export default function PredictButton() {
    const { state, dispatch } = useEditor();

    const handlePredict = async () => {
        dispatch({ type: "SET_BUSY", payload: true });
        try {
            const payload = {
                atoms: state.atoms.map((a: any) => ({ id: a.id, element: a.element, x: a.position[0], y: a.position[1], z: a.position[2] || 0 })),
                bonds: state.bonds.map((b: any) => ({ a: b.a, b: b.b, order: b.order || 1 }))
            };
            const res = await MLAPI.predictProperties(payload);
            dispatch({ type: "SET_PREDICTIONS", payload: res.properties });
        } catch (e) {
            console.warn("predict failed", e);
        } finally {
            dispatch({ type: "SET_BUSY", payload: false });
        }
    };

    return <button className="px-2 py-1 bg-blue-600 text-white rounded hover:bg-blue-700 disabled:opacity-50" onClick={handlePredict} disabled={state.busy}>Predict Properties</button>;
}
