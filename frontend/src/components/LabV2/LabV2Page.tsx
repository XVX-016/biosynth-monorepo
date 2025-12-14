import TopActionBar from "./layout/TopActionBar";
import LeftToolDock from "./layout/LeftToolDock";
import RightInspector from "./layout/RightInspector";
import BottomTemplateBar from "./layout/BottomTemplateBar";
import LabCanvas from "./LabCanvas";
import Navbar from "../Navbar";

export default function LabV2Page() {
    return (
        <div className="h-screen w-full grid grid-rows-[64px_64px_1fr_72px] overflow-hidden bg-white">
            {/* Row 1: Global Navbar */}
            <div className="z-30 border-b border-gray-100 bg-white">
                <Navbar />
            </div>

            {/* Row 2: Top Toolbar */}
            <div className="z-20 border-b border-gray-100 bg-white">
                <TopActionBar />
            </div>

            {/* Row 3: Body (Palette + Canvas) */}
            <div className="grid grid-cols-[72px_1fr] min-h-0 relative z-0">
                {/* Col 1: Palette */}
                <div className="h-full z-10">
                    <LeftToolDock />
                </div>

                {/* Col 2: Canvas */}
                <div className="relative w-full h-full bg-[#f8f9fa] overflow-hidden">
                    {/* The Sacred Canvas */}
                    <div className="absolute inset-0 z-0">
                        <LabCanvas />
                    </div>

                    {/* Floating Inspector (Overlay) */}
                    {/* We keep this here to slide over canvas without resizing it, 
                        preserving the 'Maximum uninterrupted canvas' goal until interaction */}
                    <RightInspector />
                </div>
            </div>

            {/* Row 4: Bottom Templates */}
            <div className="z-20 border-t border-gray-100 bg-white">
                <BottomTemplateBar />
            </div>
        </div>
    );
}
