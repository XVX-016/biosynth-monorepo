import { LabCanvas } from "../../components/lab/LabCanvas";
import { LabToolbar } from "../../components/lab/LabToolbar";

export default function Lab() {
    return (
        <div className="w-full h-[calc(100vh-80px)] flex bg-white">
            {/* Left Sidebar Toolbar */}
            <div className="w-20 border-r border-gray-200 bg-white">
                <LabToolbar />
            </div>

            {/* 3D Canvas - Full Area */}
            <div className="flex-1">
                <LabCanvas />
            </div>
        </div>
    );
}
