import { EditorProvider } from "../../context/EditorContext";
import { motion } from "framer-motion";
import FloatingToolbar from "../../components/LabV2/FloatingToolbar";
import LeftPanel from "../../components/LabV2/LeftPanel";
import RightPanel from "../../components/LabV2/RightPanel";
import LabCanvas from "../../components/LabV2/LabCanvas";

export default function LabV2Page() {
    return (
        <EditorProvider>
            <motion.div
                initial={{ opacity: 0 }}
                animate={{ opacity: 1 }}
                transition={{ duration: 0.3 }}
                style={{ width: "100vw", height: "100vh", position: "relative", background: "#0b0c0f" }}
            >
                <FloatingToolbar />
                <LeftPanel />
                <RightPanel />
                <LabCanvas />
            </motion.div>
        </EditorProvider>
    );
}
