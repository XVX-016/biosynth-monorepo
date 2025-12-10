import { useLabStore } from "../../store/labStore";

export default function PropertiesPanel() {
    const { selectedAtomId, selectedBondId, molecule } = useLabStore();

    const selectedAtom = selectedAtomId ? molecule.atoms.find(a => a.id === selectedAtomId) : null;
    const selectedBond = selectedBondId ? molecule.bonds.find(b => b.id === selectedBondId) : null;

    return (
        <div className="bg-white border border-gray-200 shadow-md rounded-xl p-4 w-64 pointer-events-auto">
            <h2 className="text-gray-700 font-semibold mb-3 text-sm">Properties</h2>
            <div className="space-y-3 text-gray-600 text-sm">
                {selectedAtom && (
                    <div className="space-y-1">
                        <div className="font-medium text-gray-800 border-b pb-1 mb-2">Atom</div>
                        <div className="flex justify-between"><span>Element:</span> <span className="font-mono">{selectedAtom.element}</span></div>
                        <div className="flex justify-between"><span>ID:</span> <span className="font-mono text-xs">{selectedAtom.id.slice(0, 4)}...</span></div>
                        <div className="mt-2 text-xs text-gray-500">
                            Pos: [{selectedAtom.position.map(v => v.toFixed(2)).join(', ')}]
                        </div>
                    </div>
                )}
                {selectedBond && (
                    <div className="space-y-1">
                        <div className="font-medium text-gray-800 border-b pb-1 mb-2">Bond</div>
                        <div className="flex justify-between"><span>Order:</span> <span className="font-mono">{selectedBond.order}</span></div>
                        <div className="flex justify-between"><span>Type:</span> <span className="font-mono">{selectedBond.order === 1 ? 'Single' : selectedBond.order === 2 ? 'Double' : 'Triple'}</span></div>
                    </div>
                )}
                {!selectedAtom && !selectedBond && (
                    <div className="text-gray-400 italic py-4 text-center">No atom or bond selected</div>
                )}
            </div>
        </div>
    );
}
