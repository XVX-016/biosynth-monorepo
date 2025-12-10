import type { Atom, Bond } from '../../context/EditorContext';

// small orchestration layer: call backend ML when we need bond predictions
export default class BondEngine {
    baseUrl: string;
    constructor(baseUrl: string = '/api') { this.baseUrl = baseUrl; }

    async predictBonds(atoms: Atom[]): Promise<Bond[]> {
        // send minimal JSON to backend
        const payload = { atoms: atoms.map(a => ({ id: a.id, element: a.element, position: a.position })) };
        const res = await fetch(`${this.baseUrl}/lab/predict-bonds`, {
            method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(payload)
        });
        if (!res.ok) throw new Error('Bond prediction failed');
        const json = await res.json();
        return json.bonds as Bond[];
    }
}
