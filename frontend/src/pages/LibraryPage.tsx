import MoleculeCard from "../components/Library/MoleculeCard";
import { useLibrary } from "../hooks/useLibrary";

export default function LibraryPage() {
    const { molecules, loading } = useLibrary();

    return (
        <div className="flex w-full h-[calc(100vh-64px)]">
            {/* SIDEBAR */}
            <div className="w-64 p-4 border-r bg-gray-50 flex-shrink-0">
                <input
                    className="w-full border rounded px-3 py-2 mb-4"
                    placeholder="Search molecules..."
                />

                <h3 className="font-semibold text-gray-700 mb-2">Filters</h3>
                <div className="flex flex-col gap-2 mb-6">
                    <button className="filter-btn text-left">Organic</button>
                    <button className="filter-btn text-left">Inorganic</button>
                    <button className="filter-btn text-left">Valid</button>
                </div>

                <h3 className="font-semibold text-gray-700 mb-2">Upload Molecule</h3>
                <div className="bg-white border rounded p-3">
                    <p className="text-xs text-gray-500 mb-2">Drag & drop .json files</p>
                    <button className="w-full bg-blue-500 text-white py-2 rounded text-sm hover:bg-blue-600">
                        Upload
                    </button>
                </div>
            </div>

            {/* GRID */}
            <div className="flex-1 p-6 overflow-y-auto bg-white">
                {loading ? (
                    <div className="text-gray-500">Loading library...</div>
                ) : molecules.length === 0 ? (
                    <div className="text-gray-500">No molecules found. Upload one to get started!</div>
                ) : (
                    <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
                        {molecules.map(m => (
                            <MoleculeCard key={m.id} mol={m} />
                        ))}
                    </div>
                )}
            </div>
        </div>
    );
}
