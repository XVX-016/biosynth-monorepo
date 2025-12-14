import { useLabStore } from "../../../store/labStore";

export default function RightInspector() {
    const { selectedAtomId, selectedBondId, molecule } = useLabStore();

    const selectedAtom = selectedAtomId ? molecule.atoms.find(a => a.id === selectedAtomId) : null;
    const selectedBond = selectedBondId ? molecule.bonds.find(b => b.id === selectedBondId) : null;

    // Only show if something is selected
    const isVisible = !!selectedAtom || !!selectedBond;

    return (
        <div
            className={`
                fixed right-4 top-28 bottom-24
                w-80
                bg-white/95 backdrop-blur
                rounded-2xl shadow-lg
                transition-transform duration-300 ease-in-out
                flex flex-col overflow-hidden
                z-50
                border border-gray-100
                ${isVisible ? 'translate-x-0 opacity-100' : 'translate-x-[110%] opacity-0 pointer-events-none'}
            `}
        >
            <div className="p-4 border-b border-gray-100 flex justify-between items-center bg-gray-50/50">
                <h3 className="font-semibold text-gray-800">Inspector</h3>
                <span className="text-xs text-gray-400 uppercase tracking-wider">{selectedAtom ? 'Atom' : selectedBond ? 'Bond' : 'Properties'}</span>
            </div>

            <div className="flex-1 p-4 overflow-y-auto">
                {selectedAtom && (
                    <div className="space-y-4">
                        <div className="flex items-center gap-3 mb-2">
                            <div className="w-12 h-12 rounded-xl bg-gray-100 flex items-center justify-center text-xl font-bold border border-gray-200">
                                {selectedAtom.element}
                            </div>
                            <div>
                                <div className="font-semibold text-gray-900">{selectedAtom.element} Atom</div>
                                <div className="text-xs text-gray-500 font-mono">ID: {selectedAtom.id.slice(0, 8)}...</div>
                            </div>
                        </div>

                        <div className="bg-gray-50 rounded-lg p-3 border border-gray-100">
                            <h4 className="text-xs font-semibold text-gray-400 uppercase mb-2">Position</h4>
                            <div className="grid grid-cols-3 gap-2 text-sm font-mono text-gray-600">
                                <div className="flex flex-col"><span className="text-[10px] text-gray-400">X</span>{selectedAtom.position?.x.toFixed(3) ?? 'N/A'}</div>
                                <div className="flex flex-col"><span className="text-[10px] text-gray-400">Y</span>{selectedAtom.position?.y.toFixed(3) ?? 'N/A'}</div>
                                <div className="flex flex-col"><span className="text-[10px] text-gray-400">Z</span>{selectedAtom.position?.z.toFixed(3) ?? 'N/A'}</div>
                            </div>
                        </div>
                    </div>
                )}

                {selectedBond && (
                    <div className="space-y-4">
                        <div className="flex items-center gap-3 mb-2">
                            <div className="w-12 h-12 rounded-xl bg-gray-100 flex items-center justify-center font-bold border border-gray-200">
                                {selectedBond.order === 1 ? '-' : selectedBond.order === 2 ? '=' : '≡'}
                            </div>
                            <div>
                                <div className="font-semibold text-gray-900">Bond</div>
                                <div className="text-xs text-gray-500 font-mono">Order: {selectedBond.order}</div>
                            </div>
                        </div>
                        <div className="text-sm text-gray-600">
                            Connects {selectedBond.from} to {selectedBond.to}
                        </div>
                    </div>
                )}

                {!isVisible && (
                    <div className="flex items-center justify-center h-full text-gray-400 text-sm italic">
                        Select an item to view properties
                    </div>
                )}
            </div>
        </div>
    );
}
