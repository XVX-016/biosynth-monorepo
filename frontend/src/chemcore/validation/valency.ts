import type { Molecule } from "../graph/molecule";
import { getElementSpec, type BondOrder } from "../../utils/elements";

/**
 * Calculate the sum of bond orders for an atom
 */
export function getAtomBondOrderSum(
    molecule: Molecule,
    atomId: string
): number {
    return molecule.bonds.reduce((sum, bond) => {
        if (bond.from === atomId || bond.to === atomId) {
            return sum + bond.order;
        }
        return sum;
    }, 0);
}

/**
 * Check if a bond can be added to an atom
 * 
 * Validates:
 * 1. Bond order doesn't exceed element's maxBondOrder
 * 2. Total valence doesn't exceed any allowed valence value
 */
export function canAddBond(
    molecule: Molecule,
    atomId: string,
    newBondOrder: BondOrder
): boolean {
    const atom = molecule.atoms.find(a => a.id === atomId);
    if (!atom) return false;

    const spec = getElementSpec(atom.element);
    if (!spec) {
        // Unknown element: allow for flexibility (could be strict later)
        return true;
    }

    // Check max bond order
    if (newBondOrder > spec.maxBondOrder) {
        return false;
    }

    // Check valency
    const used = getAtomBondOrderSum(molecule, atomId);
    const newTotal = used + newBondOrder;

    // Check if new total fits any allowed valence
    return spec.allowedValence.some(v => newTotal <= v);
}

export interface ValidationResult {
    valid: boolean;
    atomId?: string;
    max?: number;
    used?: number;
    reason?: string;
}

/**
 * Validate valency for all atoms in a molecule
 */
export function validateValency(molecule: Molecule): ValidationResult {
    const usage: Record<string, number> = {};

    // Calculate bond order sum for each atom
    for (const bond of molecule.bonds) {
        usage[bond.from] = (usage[bond.from] || 0) + bond.order;
        usage[bond.to] = (usage[bond.to] || 0) + bond.order;
    }

    // Validate each atom
    for (const atom of molecule.atoms) {
        const used = usage[atom.id] || 0;
        const spec = getElementSpec(atom.element);

        if (spec) {
            // Check if used valence fits any allowed valence
            const isValid = spec.allowedValence.some(v => used <= v);

            if (!isValid) {
                return {
                    valid: false,
                    atomId: atom.id,
                    max: Math.max(...spec.allowedValence),
                    used,
                    reason: `Valency exceeded for ${atom.element} (Max ${Math.max(...spec.allowedValence)}, Used ${used})`
                };
            }
        }
        // Unknown elements: skip validation (could be strict later)
    }

    return { valid: true };
}
