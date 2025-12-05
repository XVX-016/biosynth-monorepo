export function LabToolbar() {
    return (
        <div className="flex gap-4 px-6 py-3 items-center">
            <button className="px-4 py-2 rounded bg-black text-white">Select</button>
            <button className="px-4 py-2 rounded bg-gray-700 text-white">Atom</button>
            <button className="px-4 py-2 rounded bg-gray-700 text-white">Bond</button>
            <button className="px-4 py-2 rounded bg-gray-700 text-white">Drag</button>
            <button className="px-4 py-2 rounded bg-gray-700 text-white">Erase</button>
        </div>
    );
}
