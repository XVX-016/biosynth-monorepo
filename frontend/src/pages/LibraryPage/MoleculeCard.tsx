import { Molecule } from "../../types/molecule";

interface MoleculeCardProps {
    molecule: Molecule;
    onOpen: () => void;
    onDelete: () => void;
}

export default function MoleculeCard({ molecule, onOpen, onDelete }: MoleculeCardProps) {
    return (
        <div className="rounded-lg shadow hover:shadow-lg transition p-3 cursor-pointer bg-white flex flex-col">
            <div className="w-full h-32 bg-gray-100 rounded mb-2 overflow-hidden">
                {molecule.previewImage ? (
                    <img src={molecule.previewImage} className="w-full h-full object-cover" alt={molecule.name} />
                ) : (
                    <div className="w-full h-full flex items-center justify-center text-gray-400 text-xs">No Preview</div>
                )}
            </div>

            <div className="mt-2 font-medium truncate">{molecule.name}</div>

            <div className="text-xs text-gray-500 truncate">{molecule.formula || "No Formula"}</div>

            <div className="flex justify-between mt-auto pt-2">
                <button className="text-blue-600 text-sm hover:underline" onClick={(e) => { e.stopPropagation(); onOpen(); }}>View</button>
                <button className="text-red-500 text-sm hover:underline" onClick={(e) => { e.stopPropagation(); onDelete(); }}>Delete</button>
            </div>
        </div>
    );
}
