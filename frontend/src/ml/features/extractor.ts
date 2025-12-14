import type { Molecule } from '../../chemcore/graph/molecule';

const ATOMIC_WEIGHTS: Record<string, number> = {
    H: 1.008, C: 12.011, N: 14.007, O: 15.999, F: 18.998,
    P: 30.974, S: 32.06, Cl: 35.45, Br: 79.904, I: 126.90
};

export interface MoleculeFeatures {
    atomCounts: Record<string, number>;
    bondCounts: { single: number; double: number; triple: number; aromatic: number };
    molecularWeight: number;
    numAtoms: number;
    numBonds: number;
}

export function extractFeatures(molecule: Molecule): MoleculeFeatures {
    const counts: Record<string, number> = {};
    let mw = 0;

    for (const atom of molecule.atoms) {
        counts[atom.element] = (counts[atom.element] || 0) + 1;
        mw += ATOMIC_WEIGHTS[atom.element] || 0;
    }

    const bondCounts = { single: 0, double: 0, triple: 0, aromatic: 0 };
    for (const bond of molecule.bonds) {
        if (bond.order === 1) bondCounts.single++;
        else if (bond.order === 2) bondCounts.double++;
        else if (bond.order === 3) bondCounts.triple++;
        else bondCounts.aromatic++;
    }

    return {
        atomCounts: counts,
        bondCounts,
        molecularWeight: mw,
        numAtoms: molecule.atoms.length,
        numBonds: molecule.bonds.length
    };
}
