import { MLAPI } from "../../api/ml";
import type { MLRequest, OptimizeResponse } from "../../types/ml";

export async function runEnergyLoop(getState: () => any, applyFn: (resp: OptimizeResponse) => void, opts?: {
    maxRounds?: number, tolerance?: number
}) {
    const maxRounds = opts?.maxRounds ?? 6;
    const tolerance = opts?.tolerance ?? 1e-3;
    let prevScore = Infinity;

    for (let round = 0; round < maxRounds; round++) {
        const state = getState();
        const payload: MLRequest = {
            atoms: state.atoms.map((a: any) => ({ id: a.id, element: a.element, x: a.x, y: a.y, z: a.z ?? 0 })),
            bonds: state.bonds.map((b: any) => ({ a: b.a, b: b.b, order: b.order || 1 }))
        };

        let resp: OptimizeResponse;
        try {
            resp = await MLAPI.optimize(payload);
        } catch (e) {
            console.warn("ML optimize failed", e);
            return { success: false, reason: "ml-failed" };
        }

        applyFn(resp);
        if (Math.abs(prevScore - resp.score) < tolerance) {
            return { success: true, rounds: round + 1, finalScore: resp.score };
        }
        prevScore = resp.score;
    }

    return { success: true, rounds: maxRounds, finalScore: prevScore };
}
