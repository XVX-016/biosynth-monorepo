import { useEditorContext } from "../../context/EditorContext";
import { CH4, BENZENE } from "../../utils/defaultMolecules";

export default function LeftPanel() {
    const { dispatch } = useEditorContext();
    const elements = ["H", "C", "N", "O", "S", "P", "F", "Cl"];

    const insertTemplate = (mol: any) => {
        const prefix = Date.now().toString().slice(-4);
        const atoms = mol.atoms.map((a: any) => ({ ...a, id: `${prefix}_${a.id}` }));
        const bonds = mol.bonds.map((b: any) => ({ ...b, id: `${prefix}_b_${b.a}_${b.b}`, a: `${prefix}_${b.a}`, b: `${prefix}_${b.b}` }));
        dispatch({ type: "SET_MOLECULE", payload: { atoms, bonds, name: mol.name } });
    };

    return (
        <div style={{ position: "absolute", left: 12, top: 90, zIndex: 70, width: 120, background: "rgba(255,255,255,0.04)", padding: 12, borderRadius: 8 }}>
            <div style={{ marginBottom: 8, color: "#ddd", fontWeight: 700 }}>Templates</div>
            <button className="lab-btn block" onClick={() => insertTemplate(CH4)}>Methane</button>
            <button className="lab-btn block mt-2" onClick={() => insertTemplate(BENZENE)}>Benzene</button>

            <hr style={{ margin: "10px 0", borderColor: "rgba(255,255,255,0.04)" }} />

            <div style={{ color: "#ddd", fontWeight: 700, marginBottom: 6 }}>Elements</div>
            <div style={{ display: "grid", gridTemplateColumns: "repeat(2,1fr)", gap: 6 }}>
                {elements.map(el => <button key={el} className="lab-btn" onClick={() => dispatch({ type: "SET_TOOL", payload: "add-atom" })}>{el}</button>)}
            </div>
        </div>
    );
}
