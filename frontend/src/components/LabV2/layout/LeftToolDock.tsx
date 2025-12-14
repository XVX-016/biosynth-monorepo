import React, { useState } from 'react';
import { useLabStore } from "../../../store/labStore";
import { FolderOpen } from "lucide-react";

const ATOMS = [
    { symbol: 'C', color: '#2B2B2B' },
    { symbol: 'H', color: '#EDEDED', textColor: 'black' },
    { symbol: 'O', color: '#E53935' },
    { symbol: 'N', color: '#1E88E5' },
    { symbol: 'S', color: '#FBC02D', textColor: 'black' },
    { symbol: 'F', color: '#43A047' },
    { symbol: 'Cl', color: '#2E7D32' },
    { symbol: 'Br', color: '#8D6E63' },
];

const BONDS = [
    { type: 1, label: 'Single', icon: <div className="h-0.5 w-4 bg-current" /> },
    { type: 2, label: 'Double', icon: <div className="flex flex-col gap-0.5"><div className="h-0.5 w-4 bg-current" /><div className="h-0.5 w-4 bg-current" /></div> },
    { type: 3, label: 'Triple', icon: <div className="flex flex-col gap-0.5"><div className="h-0.5 w-4 bg-current" /><div className="h-0.5 w-4 bg-current" /><div className="h-0.5 w-4 bg-current" /></div> },
];

export default function LeftToolDock() {
    const { currentTool, setTool, setBondOrder, setCurrentElement } = useLabStore();
    const [selectedElement, setSelectedElement] = useState('C');
    // We will wire this to store later, for now visual state
    const [selectedBondData, setSelectedBondData] = useState(1);

    return (
        <div className="w-full h-full bg-white border-r border-gray-100 flex flex-col items-center py-4 gap-6 overflow-y-auto">

            {/* Atoms */}
            <div className="flex flex-col gap-3 w-full px-2 items-center">
                <span className="text-[10px] uppercase font-bold text-gray-400 tracking-wider">Atoms</span>
                {ATOMS.map(atom => (
                    <button
                        key={atom.symbol}
                        onClick={() => {
                            setTool('add-atom');
                            setCurrentElement(atom.symbol);
                            setSelectedElement(atom.symbol);
                        }}
                        className={`
                            w-10 h-10 rounded-full flex items-center justify-center transition-all duration-200 shadow-sm border
                            ${selectedElement === atom.symbol && currentTool === 'add-atom'
                                ? 'scale-110 ring-2 ring-blue-500 ring-offset-2 border-transparent'
                                : 'hover:scale-105 border-gray-200'
                            }
                        `}
                        style={{ backgroundColor: atom.color }}
                        title={`Add ${atom.symbol}`}
                    >
                        <span
                            className="text-xs font-bold"
                            style={{ color: atom.textColor || 'white' }}
                        >
                            {atom.symbol}
                        </span>
                    </button>
                ))}
            </div>

            <div className="w-8 h-px bg-gray-200" />

            {/* Bonds */}
            <div className="flex flex-col gap-3 w-full px-2 items-center">
                <span className="text-[10px] uppercase font-bold text-gray-400 tracking-wider">Bonds</span>
                {BONDS.map(bond => (
                    <button
                        key={bond.type}
                        onClick={() => {
                            setTool('add-bond');
                            setSelectedBondData(bond.type);
                            setBondOrder(bond.type as 1 | 2 | 3);
                        }}
                        className={`
                            w-10 h-10 rounded-lg flex items-center justify-center transition-all duration-200 border
                            ${selectedBondData === bond.type && currentTool === 'add-bond'
                                ? 'bg-black text-white border-black'
                                : 'bg-white text-gray-600 border-gray-200 hover:bg-gray-50'
                            }
                        `}
                        title={`${bond.label} Bond`}
                    >
                        {bond.icon}
                    </button>
                ))}
            </div>

            <div className="flex-1" />

            {/* Footer Actions */}
            <div className="w-full px-2 flex flex-col gap-2">
                <button className="w-10 h-10 mx-auto rounded-lg flex items-center justify-center text-gray-500 hover:bg-gray-100 hover:text-gray-900 transition-colors" title="Templates">
                    <FolderOpen size={20} />
                </button>
            </div>
        </div>
    );
}
