import React from "react";

export default function ToolbarButton({ label, icon, active, onClick }: { label: string; icon?: React.ReactNode; active?: boolean; onClick?: () => void }) {
    return (
        <button
            onClick={onClick}
            className={`
        flex items-center gap-2 px-3 py-1.5
        bg-white border rounded-lg shadow-sm
        hover:bg-gray-100 transition whitespace-nowrap
        ${active ? "border-blue-500 text-blue-600" : "border-gray-300 text-gray-600"}
      `}
        >
            {icon && <span>{icon}</span>}
            <span className="text-sm font-medium">{label}</span>
        </button>
    );
}
