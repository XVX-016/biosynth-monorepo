import { useEditorContext } from "../../context/EditorContext";

export default function RightPanel() {
    const { state } = useEditorContext();
    return (
        <div style={{ position: "absolute", right: 12, top: 90, zIndex: 70, width: 220, background: "rgba(255,255,255,0.03)", padding: 12, borderRadius: 8 }}>
            <div style={{ color: "#ddd", fontWeight: 700 }}>Selection</div>
            <div style={{ color: "#bbb", fontSize: 13 }}>
                {state.selection ? `${state.selection.type} ${state.selection.id}` : "No atom or bond selected"}
            </div>

            <hr style={{ margin: "10px 0", borderColor: "rgba(255,255,255,0.04)" }} />
            <div style={{ color: "#ddd", fontWeight: 700 }}>Info</div>
            <div style={{ color: "#bbb", fontSize: 13 }}>Atoms: {state.atoms.length}</div>
            <div style={{ color: "#bbb", fontSize: 13 }}>Bonds: {state.bonds.length}</div>
            <div style={{ color: "#bbb", fontSize: 13 }}>Tool: {state.tool}</div>
        </div>
    );
}
