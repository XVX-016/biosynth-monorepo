import ToolbarButton from "./ToolbarButton";
import { useLabStore } from "../../store/labStore";
import { LibraryAPI } from "../../api/library";

export default function FloatingToolbar() {
    const { currentTool, setTool, resetMolecule, molecule } = useLabStore();

    const handleSave = async () => {
        const name = prompt("Enter molecule name", "New Molecule");
        if (!name) return;
        try {
            await LibraryAPI.upload({
                name,
                json_graph: { atoms: molecule.atoms, bonds: molecule.bonds }
            });
            alert("Saved to Library!");
        } catch (e) {
            console.error(e);
            alert("Failed to save.");
        }
    };

    return (
        <div className="flex gap-3 bg-white/80 backdrop-blur-sm p-3 rounded-xl shadow border border-gray-200">
            <ToolbarButton
                label="Select"
                active={currentTool === 'select'}
                onClick={() => setTool('select')}
            />

            <ToolbarButton
                label="Add Atom"
                active={currentTool === 'add-atom'}
                onClick={() => setTool('add-atom')}
            />

            <ToolbarButton
                label="Add Bond"
                active={currentTool === 'add-bond'}
                onClick={() => setTool('add-bond')}
            />

            <ToolbarButton
                label="AutoBond"
                onClick={() => console.log("AutoBond triggered")}
            />

            <ToolbarButton
                label="Optimize"
                onClick={() => console.log("Optimize triggered")}
            />

            <ToolbarButton
                label="Clear"
                onClick={() => resetMolecule()}
            />

            <div className="w-px bg-gray-300 mx-1"></div>

            <ToolbarButton
                label="Save"
                onClick={handleSave}
            />
        </div>
    );
}
