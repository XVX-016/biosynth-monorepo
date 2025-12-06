import MoleculeCard from "./MoleculeCard";
import MoleculeViewerModal from "./MoleculeViewerModal";
import { useState } from "react";
import { Molecule } from "../../types/molecule";

interface LibraryContentProps {
    molecules: Molecule[];
    loading: boolean;
    remove: (id: string) => Promise<void>;
}

export default function LibraryContent({ molecules, loading, remove }: LibraryContentProps) {
    const [selected, setSelected] = useState<Molecule | null>(null);

    if (loading) return <div className="p-6">Loading...</div>;

    return (
        <div className="flex-1 p-6 grid grid-cols-4 gap-4 overflow-y-auto">
            {molecules.map(m => (
                <MoleculeCard
                    key={m.id}
                    molecule={m}
                    onOpen={() => setSelected(m)}
                    onDelete={() => remove(m.id)}
                />
            ))}

            {selected && (
                <MoleculeViewerModal molecule={selected} onClose={() => setSelected(null)} />
            )}
        </div>
    );
}
