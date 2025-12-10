import { useNavigate } from "react-router-dom";
import type { Molecule } from "../../types/molecule";

export default function MoleculeCard({ mol }: { mol: Molecule }) {
    const navigate = useNavigate();

    return (
        <div className="bg-white shadow-sm rounded-lg p-4 border hover:shadow-md transition cursor-pointer" onClick={() => navigate(`/lab?load=${mol.id}`)}>
            <div className="h-32 w-full bg-gray-100 rounded mb-3 overflow-hidden">
                {mol.previewImage ? (
                    <img src={mol.previewImage} alt={mol.name} className="w-full h-full object-cover" />
                ) : (
                    <div className="w-full h-full flex items-center justify-center text-gray-400">No Preview</div>
                )}
            </div>

            <div className="font-semibold text-gray-800 mb-1">
                {mol.name}
            </div>

            <div className="flex gap-2 text-xs mb-3">
                <span className="px-2 py-1 bg-blue-50 text-blue-600 rounded">Molecule</span>
            </div>

            <div className="flex gap-2">
                <button
                    onClick={(e) => { e.stopPropagation(); navigate(`/lab?load=${mol.id}`); }}
                    className="px-3 py-1 bg-black text-white rounded text-sm"
                >
                    View in Lab
                </button>

                <button
                    onClick={(e) => { e.stopPropagation(); console.log("fork"); }}
                    className="px-3 py-1 bg-gray-200 rounded text-sm"
                >
                    Fork
                </button>
            </div>
        </div>
    );
}
