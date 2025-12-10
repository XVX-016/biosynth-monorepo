export type Atom = { id: string; element: string; x: number; y: number; z: number };
export type Bond = { id: string; a: string; b: string; order: number };

const covalentRadii: Record<string, number> = {
    H: 0.31, C: 0.76, N: 0.71, O: 0.66, F: 0.57,
    P: 1.07, S: 1.05, Cl: 1.02, Br: 1.20, I: 1.39,
};

const maxValence: Record<string, number> = {
    H: 1, C: 4, N: 3, O: 2, F: 1,
    S: 6, P: 5, Cl: 1, Br: 1, I: 1,
};

const BOND_TOLERANCE = 0.45;
const MIN_BOND_DISTANCE = 0.3;

function dist(a: Atom, b: Atom) {
    return Math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2);
}

function calculateValence(atomId: string, bonds: Bond[]) {
    return bonds
        .filter(b => b.a === atomId || b.b === atomId)
        .reduce((acc, b) => acc + b.order, 0);
}

export class BondEngine {
    static shouldBond(a: Atom, b: Atom): boolean {
        const r1 = covalentRadii[a.element];
        const r2 = covalentRadii[b.element];
        if (!r1 || !r2) return false;
        const ideal = r1 + r2;
        const d = dist(a, b);
        return d > MIN_BOND_DISTANCE && d < ideal + BOND_TOLERANCE && d < ideal * 1.3;
    }

    static canAddBond(a: Atom, b: Atom, bonds: Bond[]) {
        const maxA = maxValence[a.element] ?? 0;
        const maxB = maxValence[b.element] ?? 0;
        const valA = calculateValence(a.id, bonds);
        const valB = calculateValence(b.id, bonds);
        if (valA >= maxA || valB >= maxB) return false;
        if (bonds.some(bd => (bd.a === a.id && bd.b === b.id) || (bd.a === b.id && bd.b === a.id))) return false;
        return true;
    }

    static autoBond(atoms: Atom[], bonds: Bond[]) {
        const newBonds: Bond[] = [];
        const temp = bonds.slice();
        for (let i = 0; i < atoms.length; i++) {
            for (let j = i + 1; j < atoms.length; j++) {
                const a = atoms[i];
                const b = atoms[j];
                if (this.shouldBond(a, b) && this.canAddBond(a, b, temp)) {
                    const nb: Bond = { id: crypto.randomUUID(), a: a.id, b: b.id, order: 1 };
                    newBonds.push(nb);
                    temp.push(nb);
                }
            }
        }
        return newBonds;
    }
}
