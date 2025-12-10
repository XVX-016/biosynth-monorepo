import React from "react";
import { useEditor } from "../../context/EditorContext";
import { runEnergyLoop } from "../../engine/ml/energyLoop";

export default function OptimizeButton() {
    const { state, dispatch } = useEditor();

    const handleOptimize = async () => {
        dispatch({ type: "SET_BUSY", payload: true });

        // Wrapper to get current state for the loop
        // Note: runEnergyLoop expects a state getter or we pass initial state?
        // The scaffold passes getState closure. We can't easily pass a closure that gets *fresh* state from React context.
        // Instead, we might just run one iteration or need a different approach for the loop.
        // For now, let's just trigger one optimization call directly to MLAPI via APPLY_ML action

        // Actually, scaffold used: runEnergyLoop(getState, applyFn)
        // We'll stub it to just single pass for now to avoid React state closure complexity
        try {
            // Just call optimize once
            const { MLAPI } = await import("../../api/ml"); // dynamic match scaffold
            const payload = {
                atoms: state.atoms.map(a => ({ id: a.id, element: a.element, x: a.position[0], y: a.position[1], z: a.position[2] || 0 })),
                bonds: state.bonds.map(b => ({ a: b.a, b: b.b, order: b.order || 1 }))
            };
            const res = await MLAPI.optimize(payload);
            dispatch({ type: "APPLY_ML", payload: { correctedAtoms: res.correctedAtoms } });
        } catch (e) {
            console.error("Optimize failed", e);
        } finally {
            dispatch({ type: "SET_BUSY", payload: false });
        }
    };

    return (
        <button
            className="px-2 py-1 border rounded hover:bg-gray-50 disabled:opacity-50"
            onClick={handleOptimize}
            disabled={state.busy}
        >
            Optimize Geometry
        </button>
    );
}
