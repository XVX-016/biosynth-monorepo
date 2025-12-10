import React from 'react';
import Lab3DViewport from './Lab3DViewport';
import LabSidebar from "./LabSidebar";
import LabToolbar from "./LabToolbar";
import LabInspector from "./LabInspector";
import LabBottomPanel from "./LabBottomPanel";

export default function LabLayout() {
    return (
        <div className="flex h-screen w-full overflow-hidden">
            {/* Left Sidebar */}
            <LabSidebar />

            {/* Center Main Area */}
            <div className="flex flex-1 flex-col relative min-w-0">
                {/* Toolbar */}
                <LabToolbar />

                {/* Viewport */}
                <div className="flex-1 relative">
                    <Lab3DViewport />
                </div>

                {/* Bottom Panel */}
                <div className="h-[160px] flex-shrink-0">
                    <LabBottomPanel />
                </div>
            </div>

            {/* Right Inspector */}
            <LabInspector />
        </div>
    );
}
