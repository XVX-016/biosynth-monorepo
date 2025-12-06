// src/components/lab/OptimizeButton.tsx
import React from "react";
// We use moleculeStore instead of EditorContext because the app uses moleculeStore for Lab.
// We must bridge the gap or migrate. The prompt asked for EditorContext but the app uses Store.
// I will use useMoleculeStore for state, and adapt the loop helper.
// Actually, optimizing updates positions. 
// useMoleculeStore has updatePosition but for single atoms.
// It has setMolecule.
// I will implement optimization using store.

import { useMoleculeStore } from "../../store/moleculeStore";
import { runEnergyLoop } from "../../engine/ml/energyLoop";
import type { OptimizeResponse } from "../../types/ml";
import { MoleculeGraph } from "@biosynth/engine";

export default function OptimizeButton() {
    const currentMolecule = useMoleculeStore(state => state.currentMolecule);
    const setMolecule = useMoleculeStore(state => state.setMolecule);
    const optimizationStatus = useMoleculeStore(state => state.optimizationStatus);
    // Reuse existing status or add new logic

    const [busy, setBusy] = React.useState(false);

    const applyResp = (resp: OptimizeResponse) => {
        if (!currentMolecule) return;

        // Clone molecule and update positions
        const cloned = currentMolecule.clone();
        resp.correctedAtoms.forEach(correction => {
            const atom = cloned.atoms.get(correction.id);
            if (atom) {
                cloned.atoms.set(correction.id, {
                    ...atom,
                    position: [correction.x, correction.y, correction.z]
                });
            }
        });
        setMolecule(cloned);
    };

    const handleOptimize = async () => {
        if (!currentMolecule) return;
        setBusy(true);

        try {
            await runEnergyLoop(
                () => {
                    // Convert graph to atom/bond arrays for ML request
                    const atoms = Array.from(currentMolecule.atoms.values()).map(a => ({
                        id: a.id,
                        element: a.element,
                        x: a.position[0],
                        y: a.position[1],
                        z: a.position[2]
                    }));
                    const bonds = Array.from(currentMolecule.bonds.values()).map(b => ({
                        a: b.a1,
                        b: b.a2,
                        order: b.order
                    }));
                    return { atoms, bonds };
                },
                (resp) => applyResp(resp),
                { maxRounds: 4 }
            );
        } catch (e) {
            console.error("Optimization error", e);
        } finally {
            setBusy(false);
        }
    };

    return (
        <button
            onClick={handleOptimize}
            disabled={busy || !currentMolecule}
            className={`px-3 py-1 bg-blue-600 text-white rounded text-sm hover:bg-blue-700 disabled:opacity-50 ${busy ? "animate-pulse" : ""}`}
        >
            {busy ? "Optimizing..." : "Optimize Geometry"}
        </button>
    );
}
