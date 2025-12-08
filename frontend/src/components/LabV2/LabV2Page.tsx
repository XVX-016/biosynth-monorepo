import FloatingToolbar from "./FloatingToolbar";
import LabCanvas from "./LabCanvas";
import TemplatesPanel from "./TemplatesPanel";
import { EditorProvider } from "../../context/EditorContext";

/**
 * LabV2Page
 * Wraps EditorProvider and renders floating toolbar + canvas + right inspector.
 */
export default function LabV2Page() {
    return (
        <EditorProvider>
            <div className="w-full h-screen relative bg-gray-900">
                {/* Floating toolbar (draggable) */}
                <FloatingToolbar />

                {/* Left templates / quick actions (small translucent panel) */}
                <div style={{
                    position: "absolute",
                    left: 16,
                    top: 96,
                    zIndex: 60,
                }}>
                    <TemplatesPanel />
                </div>

                {/* Main 3D canvas */}
                <div style={{ width: "100%", height: "100%" }}>
                    <LabCanvas />
                </div>

                {/* Right inspector */}
                <div style={{
                    position: "absolute",
                    right: 16,
                    top: 96,
                    zIndex: 60,
                }}>
                    <div className="lab-panel w-64">
                        <div className="font-semibold mb-2">Selection</div>
                        <div className="text-sm text-gray-600">No atom or bond selected</div>
                        <hr className="my-3" />
                        <div className="font-semibold">Atom Properties</div>
                        <div className="text-sm text-gray-600">Select an atom</div>
                        <hr className="my-3" />
                        <div className="font-semibold">Bond Properties</div>
                        <div className="text-sm text-gray-600">Select a bond</div>
                    </div>
                </div>

            </div>
        </EditorProvider>
    );
}
