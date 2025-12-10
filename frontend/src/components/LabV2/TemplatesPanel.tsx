import { CH4, BENZENE } from "../../utils/defaultMolecules";
import { useLabStore } from "../../store/labStore";
import { nanoid } from "nanoid";

export default function TemplatesPanel() {
    const { loadMolecule } = useLabStore();

    const insert = (template: any) => {
        // Map template (x,y,z fields) to Store (position array)
        // And remap IDs to avoid collisions
        const idMap: Record<string, string> = {};

        const atoms = template.atoms.map((a: any) => {
            const newId = nanoid();
            idMap[a.id] = newId;
            return {
                id: newId,
                element: a.element,
                position: [a.x, a.y, a.z || 0],
                charge: 0
            };
        });

        const bonds = template.bonds.map((b: any) => ({
            id: nanoid(),
            atom1: idMap[b.a],
            atom2: idMap[b.b],
            order: b.order || 1
        }));

        loadMolecule({
            id: nanoid(),
            name: template.name,
            atoms,
            bonds,
            formula: "",
            createdAt: new Date().toISOString(),
            updatedAt: new Date().toISOString(),
            isValid: true,
            qualityScore: 100,
            previewImage: ""
        });
    };

    return (
        <div className="bg-white border border-gray-200 shadow-md rounded-xl p-4 w-64 pointer-events-auto">
            <h2 className="text-gray-700 font-semibold mb-3 text-sm">Templates</h2>
            <div className="space-y-2 text-gray-600 text-sm">
                <button
                    onClick={() => insert(CH4)}
                    className="w-full text-left px-3 py-2 bg-gray-50 hover:bg-gray-100 rounded text-gray-700 transition"
                >
                    Methane
                </button>
                <button
                    onClick={() => insert(BENZENE)}
                    className="w-full text-left px-3 py-2 bg-gray-50 hover:bg-gray-100 rounded text-gray-700 transition"
                >
                    Benzene
                </button>
            </div>
        </div>
    );
}
