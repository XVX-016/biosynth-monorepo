import { useNavigate } from "react-router-dom";
import type { Molecule } from "../../types/molecule";

interface MoleculeCardProps {
    mol: Molecule;
    onView?: () => void;
    onFork?: () => void;
}

export default function MoleculeCard({ mol, onView, onFork }: MoleculeCardProps) {
    const navigate = useNavigate();

    return (
        <div
            className="bg-white shadow-sm rounded-lg p-4 border hover:shadow-md transition cursor-pointer group flex flex-col"
            onClick={onView}
        >
            <div className="h-32 w-full bg-gray-100 rounded mb-3 overflow-hidden relative">
                {mol.previewImage ? (
                    <img src={mol.previewImage} alt={mol.name} className="w-full h-full object-cover" />
                ) : (
                    <div className="w-full h-full flex items-center justify-center text-gray-400">
                        <span className="text-4xl text-gray-200">⬡</span>
                    </div>
                )}

                {/* Quick actions on hover */}
                <div className="absolute inset-0 bg-black/10 opacity-0 group-hover:opacity-100 transition flex items-center justify-center gap-2">
                    <button
                        onClick={(e) => { e.stopPropagation(); navigate(`/lab?load=${mol.id}`); }}
                        className="bg-white text-black text-xs px-3 py-1.5 rounded-full font-medium hover:scale-105 transition shadow-sm"
                    >
                        Open in Lab
                    </button>
                </div>
            </div>

            <div className="font-semibold text-gray-800 mb-1 truncate" title={mol.name}>
                {mol.name}
            </div>

            <div className="flex gap-2 text-xs mb-3">
                <span className="px-2 py-0.5 bg-blue-50 text-blue-600 rounded border border-blue-100">Molecule</span>
                {mol.isValid && <span className="px-2 py-0.5 bg-green-50 text-green-600 rounded border border-green-100">Valid</span>}
            </div>

            <div className="flex gap-2 mt-auto">
                <button
                    onClick={(e) => { e.stopPropagation(); onView?.(); }}
                    className="flex-1 px-3 py-1.5 bg-white border border-gray-200 hover:bg-gray-50 text-gray-700 rounded text-sm transition"
                >
                    View Details
                </button>

                <button
                    onClick={(e) => { e.stopPropagation(); onFork?.(); }}
                    className="px-3 py-1.5 bg-white border border-gray-200 hover:bg-gray-50 text-gray-700 rounded text-sm transition"
                    title="Fork to My Library"
                >
                    Fork
                </button>
            </div>
        </div>
    );
}
