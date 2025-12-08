import { CH4, BENZENE } from "../../utils/defaultMolecules";
import { useEditor } from "../../context/EditorContext";

/** Small translucent card with template buttons */
export default function TemplatesPanel() {
    const { dispatch } = useEditor();

    const insert = (mol: any) => {
        const prefix = Date.now().toString().slice(-4);
        const atoms = mol.atoms.map((a: any) => ({
            ...a,
            id: `${prefix}_${a.id}`,
            position: [a.x, a.y, a.z || 0]
        }));
        const bonds = mol.bonds.map((b: any) => ({
            ...b,
            a: `${prefix}_${b.a}`,
            b: `${prefix}_${b.b}`,
            id: `${prefix}_b_${b.a}_${b.b}`
        }));
        dispatch({ type: "LOAD_MOLECULE", payload: { atoms, bonds, name: mol.name } });
    };

    return (
        <div className="lab-panel">
            <div className="font-semibold mb-2">Templates</div>
            <button className="lab-btn block w-full mb-2" onClick={() => insert(CH4)}>
                Methane
            </button>
            <button className="lab-btn block w-full" onClick={() => insert(BENZENE)}>
                Benzene
            </button>
        </div>
    );
}
