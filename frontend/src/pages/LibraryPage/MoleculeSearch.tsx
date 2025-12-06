import { Search } from "lucide-react";

export default function MoleculeSearch() {
    return (
        <div className="relative">
            <Search className="absolute left-3 top-1/2 -translate-y-1/2 text-gray-400" size={16} />
            <input
                className="w-full pl-9 pr-3 py-2 rounded border border-gray-300 text-sm focus:outline-none focus:ring-2 focus:ring-blue-500"
                placeholder="Search molecules..."
            />
        </div>
    );
}
