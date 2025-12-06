export default function MoleculeFilters() {
    return (
        <div className="border border-gray-200 rounded-lg p-3 bg-white">
            <div className="font-medium text-sm text-gray-700 mb-2">Filters</div>
            <div className="flex flex-wrap gap-2">
                <button className="px-3 py-1 bg-gray-100 hover:bg-gray-200 rounded-full text-xs font-medium text-gray-600">Organic</button>
                <button className="px-3 py-1 bg-gray-100 hover:bg-gray-200 rounded-full text-xs font-medium text-gray-600">Inorganic</button>
                <button className="px-3 py-1 bg-gray-100 hover:bg-gray-200 rounded-full text-xs font-medium text-gray-600">Valid</button>
            </div>
        </div>
    );
}
