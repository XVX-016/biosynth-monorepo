import { useEffect, useState } from "react";
import { LibraryAPI } from "../api/library";
import type { Molecule } from "../types/molecule";

export function useLibrary() {
    const [molecules, setMolecules] = useState<Molecule[]>([]);
    const [loading, setLoading] = useState(true);
    const [page, setPage] = useState(1);

    const load = async (p: number = 1) => {
        setLoading(true);
        try {
            const data = await LibraryAPI.list({ page: p, limit: 9 });
            setMolecules(Array.isArray(data) ? data : []);
            setPage(p);
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
        load(1);
    }, []);

    return { molecules, loading, load, upload, remove, page };
}
