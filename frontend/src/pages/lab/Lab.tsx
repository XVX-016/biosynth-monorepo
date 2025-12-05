import { LabCanvas } from "../../components/lab/LabCanvas";
import { LabToolbar } from "../../components/lab/LabToolbar";

export default function Lab() {
    return (
        <div className="w-full h-[calc(100vh-80px)] flex flex-col bg-white">
            {/* Toolbar */}
            <div className="w-full border-b bg-white shadow-sm">
                <LabToolbar />
            </div>

            {/* 3D Editor */}
            <div className="flex-1">
                <LabCanvas />
            </div>
        </div>
    );
}
