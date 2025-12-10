import ToolbarButton from "./ToolbarButton";
import { useLabStore } from "../../store/labStore";

export default function FloatingToolbar() {
    const { currentTool, setTool, resetMolecule } = useLabStore();

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
        </div>
    );
}
