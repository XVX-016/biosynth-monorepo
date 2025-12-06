import React from "react";
import { useEditorContext } from "../../context/EditorContext";
import { MLAPI } from "../../api/ml";

export default function PredictButton() {
    const { state, dispatch } = useEditorContext();

    const handlePredict = async () => {
        dispatch({ type: "SET_BUSY", payload: true });
        try {
            const payload = {
                atoms: state.atoms.map((a: any) => ({ id: a.id, element: a.element, x: a.x, y: a.y, z: a.z || 0 })),
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

    return <button className="btn" onClick={handlePredict}>Predict properties</button>;
}
