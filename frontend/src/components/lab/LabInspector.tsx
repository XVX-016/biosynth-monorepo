import React from "react";

export default function LabInspector() {
    return (
        <div className="w-[260px] border-l border-gray-300 bg-white h-full p-4 flex flex-col gap-6">

            <div>
                <h2 className="font-semibold">Selection</h2>
                <p className="text-sm text-gray-600">No atom or bond selected</p>
            </div>

            <div>
                <h2 className="font-semibold">Atom Properties</h2>
                <p className="text-sm text-gray-600 italic">Select an atom</p>
            </div>

            <div>
                <h2 className="font-semibold">Bond Properties</h2>
                <p className="text-sm text-gray-600 italic">Select a bond</p>
            </div>

        </div>
    );
}
