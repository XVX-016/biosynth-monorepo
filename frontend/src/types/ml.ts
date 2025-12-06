export type Atom = { id: string; element: string; x: number; y: number; z?: number };
export type Bond = { a: string; b: string; order: number };

export type MLRequest = { atoms: Atom[]; bonds: Bond[] };

export type OptimizeResponse = {
    correctedAtoms: { id: string; x: number; y: number; z: number }[];
    score: number;
    iterations: number;
};

export type BondOrderResponse = {
    bondOrders: { a: string; b: string; order: number; confidence: number }[];
};

export type AutoBondResponse = {
    bonds: { a: string; b: string; order: number }[];
};
