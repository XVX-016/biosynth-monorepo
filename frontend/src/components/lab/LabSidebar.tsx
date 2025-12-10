import React from "react";

export default function LabSidebar() {
    return (
        <div className="w-[240px] border-r border-gray-300 bg-white h-full p-4 flex flex-col gap-6">

            <div>
                <h2 className="font-semibold mb-2">Library</h2>
                <button className="px-2 py-1 rounded bg-black text-white text-sm w-full">
                    Open Library
                </button>
            </div>

            <div>
                <h2 className="font-semibold mb-2">Tools</h2>
                <div className="flex flex-col gap-2 text-sm">
                    <button className="border p-2 rounded hover:bg-gray-50">Select</button>
                    <button className="border p-2 rounded hover:bg-gray-50">Add Atom</button>
                    <button className="border p-2 rounded hover:bg-gray-50">Add Bond</button>
                </div>
            </div>

            <div>
                <h2 className="font-semibold mb-2">Atoms</h2>
                <div className="grid grid-cols-3 gap-2 text-center text-sm">
                    {["H", "C", "N", "O", "S", "P", "F", "Cl"].map((el) => (
                        <button
                            key={el}
                            className="border p-2 rounded hover:bg-gray-100"
                        >
                            {el}
                        </button>
                    ))}
                </div>
            </div>

        </div>
    );
}
