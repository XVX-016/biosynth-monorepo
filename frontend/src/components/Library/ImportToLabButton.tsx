import React from "react";
// Assuming api is exported properly from api/library.ts or axios.ts
// User used "import api from ../../api/axios" but our api is LibraryAPI exported from api/library.ts
// We should check what src/api/axios.ts is or use existing axios imports.
// I will create src/api/axios.ts if missing or use generic axios.
import axios from "axios";
import { useEditorContext } from "../../context/EditorContext";
import { useNavigate } from "react-router-dom";

// Helper if api module not present as user assumed
const api = axios.create({ baseURL: "/api" });

export default function ImportToLabButton({ moleculeId }: { moleculeId: string }) {
    const { dispatch } = useEditorContext();
    const navigate = useNavigate();

    const handleImport = async () => {
        try {
            // route is /molecules/{id}/export since LibraryAPI defined base as /molecules in api/library.ts
            // But user scaffold uses api.post(...) with /library prefix?
            // In src/api/library.ts we saw LibraryAPI calls `/molecules`.
            // Backend mounts it at `/api/molecules` or similar.
            // If user provided backend code mounts `library` router, let's verify path.
            // User said: "backend router apps/backend/app/routers/library.py already supports /export".
            // And in step 228 app.py: `app.include_router(library_router.router, tags=["molecules"])`.
            // library.py has NO prefix defined in step 175.
            // So it is mounted at root `/` + endpoint path? OR `/molecules`?
            // In library.py step 175: `router = APIRouter(prefix="/molecules")`. Wait, step 175 diff showed `router = APIRouter(prefix="/molecules")` was NOT in the diff, but earlier view showed it?
            // Step 156 view of library.py showed `router = APIRouter(prefix="/molecules")`.
            // So path is `/molecules/{id}/export`.

            const res = await api.post(`/molecules/${moleculeId}/export`);
            const mol = res.data.data ?? res.data; // match backend shape
            dispatch({ type: "SET_MOLECULE", payload: mol });
            // fit camera will run in Lab viewer effect
            navigate("/lab");
        } catch (e) {
            console.error("Import failed", e);
        }
    };

    return <button className="btn" onClick={handleImport}>Open in Lab</button>;
}
