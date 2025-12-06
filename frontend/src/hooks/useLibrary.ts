import { useEffect, useState } from "react";
import { LibraryAPI } from "../api/library";
import { Molecule } from "../types/molecule";

export function useLibrary() {
    const [molecules, setMolecules] = useState<Molecule[]>([]);
    const [loading, setLoading] = useState(true);

    const load = async () => {
        setLoading(true);
        try {
            const data = await LibraryAPI.list();
            setMolecules(Array.isArray(data) ? data : []);
        } catch (e) {
            console.error("Failed to load library", e);
            setMolecules([]);
        } finally {
            setLoading(false);
        }
    };

    const remove = async (id: string) => {
        try {
            await LibraryAPI.delete(id);
            setMolecules(molecules.filter(m => m.id !== id));
        } catch (e) {
            console.error("Failed to delete molecule", e);
        }
    };

    const upload = async (mol: any) => {
        try {
            const created = await LibraryAPI.upload(mol);
            setMolecules([created, ...molecules]);
        } catch (e) {
            console.error("Failed to upload molecule", e);
        }
    };

    useEffect(() => {
        load();
    }, []);

    return { molecules, loading, load, upload, remove };
}
