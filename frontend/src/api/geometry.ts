import axios from 'axios';

// Define the shape of the molecule data expected by the backend
interface MoleculeData {
    atoms: any[];
    bonds: any[];
}

export async function optimizeGeometry(molecule: MoleculeData) {
    try {
        const { data } = await axios.post("/api/ml/optimize", molecule);
        return data; // returns new atom coordinates
    } catch (error) {
        console.error("Geometry optimization failed:", error);
        throw error;
    }
}

export async function validateBond(molecule: MoleculeData, bond: any) {
    try {
        const { data } = await axios.post("/api/ml/validate-bond", { molecule, bond });
        return data;
    } catch (error) {
        console.error("Bond validation failed:", error);
        throw error;
    }
}
