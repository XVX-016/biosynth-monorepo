import React from 'react';
import { useMoleculeStore } from '../../store/moleculeStore';
import { CH4, BENZENE, WATER, AMMONIA } from '../../utils/defaultMolecules';
import { addAtom } from '../../lib/engineAdapter'; // Or manual molecule construction

export default function TemplatesPanel() {
    const setMolecule = useMoleculeStore((state) => state.setMolecule);
    const setTool = useMoleculeStore((state) => state.setTool);

    const loadTemplate = (template: any) => {
        // Convert template JSON to MoleculeGraph
        // Assuming we rely on engineAdapter or store to parse
        // Since we don't have a direct parser here, we'll manually construct or use a helper
        // For simplicity, let's assume we can trigger a load or we have to build it

        // Quick hack: Clear and Rebuild
        // Ideally: useMoleculeStore.getState().loadJSON(template)
        // But we don't have that.

        // Let's rely on setMolecule if we can construct a MoleculeGraph
        // Or specific actions.

        // Since I can't easily import MoleculeGraph here without potentially issues, 
        // I will dispatch a custom event or use available store methods if possible.
        // Actually, let's just use the store's setMolecule if it accepts a graph.
        // I'll need to reconstruct the graph.
        // For now, let's just log or try to use a utility if available.

        // User wants "Exact commands". 
        // I'll implement a simple reconstructor here.

        import('@biosynth/engine').then(({ MoleculeGraph }) => {
            const graph = new MoleculeGraph();
            template.atoms.forEach((a: any) => {
                graph.addAtom({
                    id: a.id,
                    element: a.element,
                    position: a.position
                });
            });
            template.bonds.forEach((b: any) => {
                graph.addBond({
                    id: b.id,
                    a1: b.a,
                    a2: b.b,
                    order: b.order
                });
            });
            setMolecule(graph);
            setTool('select');
        }).catch(err => {
            console.error("Failed to load engine", err);
        });
    };

    return (
        <div className="flex flex-col gap-2 p-2">
            <div className="text-xs font-bold text-gray-500 uppercase tracking-wider mb-1">Templates</div>
            <button
                onClick={() => loadTemplate(CH4)}
                className="p-2 text-sm text-left hover:bg-gray-100 rounded flex items-center gap-2"
            >
                <div className="w-8 h-8 rounded bg-gray-200 flex items-center justify-center text-xs">CH4</div>
                <span>Methane</span>
            </button>
            <button
                onClick={() => loadTemplate(BENZENE)}
                className="p-2 text-sm text-left hover:bg-gray-100 rounded flex items-center gap-2"
            >
                <div className="w-8 h-8 rounded bg-gray-200 flex items-center justify-center text-xs">C6H6</div>
                <span>Benzene</span>
            </button>
            <button
                onClick={() => loadTemplate(WATER)}
                className="p-2 text-sm text-left hover:bg-gray-100 rounded flex items-center gap-2"
            >
                <div className="w-8 h-8 rounded bg-gray-200 flex items-center justify-center text-xs">H2O</div>
                <span>Water</span>
            </button>
            <button
                onClick={() => loadTemplate(AMMONIA)}
                className="p-2 text-sm text-left hover:bg-gray-100 rounded flex items-center gap-2"
            >
                <div className="w-8 h-8 rounded bg-gray-200 flex items-center justify-center text-xs">NH3</div>
                <span>Ammonia</span>
            </button>
        </div>
    );
}
