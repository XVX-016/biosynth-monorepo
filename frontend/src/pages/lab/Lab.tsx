import { useEffect } from "react";
import { useSearchParams } from "react-router-dom";
import { LabCanvas } from "../../components/lab/LabCanvas";
import { LabToolbar } from "../../components/lab/LabToolbar";
import { useMoleculeStore } from "../../store/moleculeStore";
import { LibraryAPI } from "../../api/library";
import { moleculeFromJSON } from "../../lib/engineAdapter";
import OptimizeButton from "../../components/lab/OptimizeButton";

export default function Lab() {
    const [searchParams] = useSearchParams();
    const loadId = searchParams.get("load") || searchParams.get("id");
    const setMolecule = useMoleculeStore(state => state.setMolecule);
    const tool = useMoleculeStore(state => state.tool);

    useEffect(() => {
        const load = async () => {
            if (loadId) {
                try {
                    // Use export endpoint for clean Lab-ready data
                    const response = await LibraryAPI.export(loadId);

                    if (response && response.data) {
                        // response.data contains the json_graph
                        // If response.data is string, parse it.
                        // But API returns dict from python, accessible as object.

                        const jsonStr = typeof response.data === 'string' ? response.data : JSON.stringify(response.data);

                        const graph = moleculeFromJSON(jsonStr);
                        if (graph) {
                            setMolecule(graph);
                            console.log("Loaded molecule:", response.name);
                        } else {
                            console.warn("Failed to create graph from data:", response.data);
                        }
                    } else if (response) {
                        // Fallback check if full response is the data?
                        const graph = moleculeFromJSON(JSON.stringify(response));
                        if (graph) {
                            setMolecule(graph);
                        }
                    }
                } catch (e) {
                    console.error("Failed to load molecule", e);
                }
            }
        };

        // Only run on mount or when ID changes
        load();
    }, [loadId, setMolecule]);

    // Cleanup on unmount or just let store persist? 
    // Usually fine to persist.

    return (
        <div className="w-full h-[calc(100vh-80px)] flex bg-white">
            {/* Left Sidebar Toolbar */}
            <div className="w-20 border-r border-gray-200 bg-white z-10 shadow-sm flex flex-col items-center pb-4">
                <LabToolbar />
                <div className="mt-auto pt-2 border-t w-full flex justify-center">
                    <OptimizeButton />
                </div>
            </div>

            {/* 3D Canvas - Full Area */}
            <div className="flex-1 relative bg-gray-50">
                <LabCanvas />

                {/* Tooltip or status overlay if needed */}
                <div className="absolute top-4 left-4 pointer-events-none text-gray-400 text-xs">
                    Tool: {tool}
                </div>
            </div>
        </div>
    );
}
