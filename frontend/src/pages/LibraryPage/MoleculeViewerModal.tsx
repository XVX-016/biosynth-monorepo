import { useEffect, useState } from "react";
import Modal from "../../components/Modal";
import MoleculeCanvas from "../../components/MoleculeCanvas";
import type { Molecule } from "../../types/molecule";
import { LibraryAPI } from "../../api/library";

interface MoleculeViewerModalProps {
    molecule: Molecule;
    onClose: () => void;
}

export default function MoleculeViewerModal({ molecule, onClose }: MoleculeViewerModalProps) {
    const [details, setDetails] = useState<any>(null);
    const [loading, setLoading] = useState(true);

    useEffect(() => {
        const fetchDetails = async () => {
            try {
                const data = await LibraryAPI.get(molecule.id);
                setDetails(data);
            } catch (e) {
                console.error("Failed to load details", e);
            } finally {
                setLoading(false);
            }
        };
        fetchDetails();
    }, [molecule.id]);

    const displayMol = details || molecule;
    // Extract atoms/bonds from json_graph if available, otherwise fallback
    const atoms = displayMol.json_graph?.atoms || displayMol.atoms || [];
    const bonds = displayMol.json_graph?.bonds || displayMol.bonds || [];

    return (
        <Modal onClose={onClose}>
            <div className="flex flex-col h-full w-full">
                <div className="p-4 border-b flex justify-between items-center">
                    <h2 className="text-xl font-bold">{molecule.name}</h2>
                    <div className="text-sm text-gray-500">{molecule.formula}</div>
                </div>

                <div className="flex-1 w-full relative bg-gray-50 min-h-[500px]">
                    {loading ? (
                        <div className="absolute inset-0 flex items-center justify-center text-gray-400">
                            Loading structure...
                        </div>
                    ) : (
                        <MoleculeCanvas atoms={atoms} bonds={bonds} />
                    )}
                </div>

                <div className="p-4 border-t bg-gray-50 flex justify-end gap-2">
                    <button className="px-4 py-2 rounded bg-gray-200 hover:bg-gray-300" onClick={onClose}>Close</button>
                    <button className="px-4 py-2 rounded bg-blue-600 text-white hover:bg-blue-700" onClick={() => window.location.href = `/lab?load=${molecule.id}`}>Open in Lab</button>
                </div>
            </div>
        </Modal>
    );
}

