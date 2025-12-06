import { MLAPI } from "../../api/ml";
import type { MLRequest } from "../../types/ml";

export async function predictBondOrders(state: { atoms: any[], bonds: any[] }) {
    const payload: MLRequest = {
        atoms: state.atoms.map(a => ({ id: a.id, element: a.element, x: a.x, y: a.y, z: a.z ?? 0 })),
        bonds: state.bonds.map(b => ({ a: b.a, b: b.b, order: b.order || 1 }))
    };

    try {
        const res = await MLAPI.bondOrder(payload);
        const out: Record<string, number> = {};
        res.bondOrders.forEach((bo: any) => {
            const keyA = `${bo.a}|${bo.b}`;
            const keyB = `${bo.b}|${bo.a}`;
            out[keyA] = bo.order;
            out[keyB] = bo.order;
        });
        return out;
    } catch (e) {
        console.warn("bondOrder failed", e);
        return {};
    }
}
