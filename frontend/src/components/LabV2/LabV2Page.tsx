import Navbar from "../Navbar";
import FloatingToolbar from "./FloatingToolbar";
import TemplatesPanel from "./TemplatesPanel";
import PropertiesPanel from "./PropertiesPanel";
import LabCanvas from "./LabCanvas";

export default function LabV2Page() {
    return (
        <div className="lab-grid-bg w-full h-screen flex flex-col overflow-hidden">
            {/* Top Navbar */}
            <div className="z-50 relative bg-white shadow-sm">
                <Navbar />
            </div>

            {/* Main Content Area */}
            <div className="flex-1 relative w-full h-full">

                {/* 3D Canvas Layer */}
                <div id="lab-canvas" className="absolute inset-0 w-full h-full z-0">
                    <LabCanvas />
                </div>

                {/* Floating Toolbar (top-center) */}
                <div className="absolute top-6 left-1/2 -translate-x-1/2 z-40 pointer-events-auto">
                    <FloatingToolbar />
                </div>

                {/* Left Templates Panel */}
                <div className="absolute top-6 left-4 z-40 pointer-events-auto">
                    <TemplatesPanel />
                </div>

                {/* Right Properties Panel */}
                <div className="absolute top-6 right-4 z-40 pointer-events-auto">
                    <PropertiesPanel />
                </div>

                {/* Event Blocker for background to let events pass to canvas */}
                <div className="absolute inset-0 pointer-events-none z-10"></div>
            </div>
        </div>
    );
}
