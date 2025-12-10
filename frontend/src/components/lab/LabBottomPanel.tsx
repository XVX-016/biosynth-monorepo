import React from "react";

export default function LabBottomPanel() {
    return (
        <div className="w-full border-t border-gray-300 bg-white h-[160px] p-2 text-sm">
            <p className="text-gray-600 font-semibold mb-2">Event Log</p>
            <div className="h-[120px] overflow-auto text-xs text-gray-700 border p-2 rounded bg-gray-50">
                {/* future console logs */}
                <p>No events yet...</p>
            </div>
        </div>
    );
}
