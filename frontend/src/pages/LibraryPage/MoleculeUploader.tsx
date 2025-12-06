import { useState } from "react";
import { Upload } from "lucide-react";

interface MoleculeUploaderProps {
    onUpload: (mol: any) => Promise<void>;
}

export default function MoleculeUploader({ onUpload }: MoleculeUploaderProps) {
    const [name, setName] = useState("");
    const [file, setFile] = useState<File | null>(null);
    const [loading, setLoading] = useState(false);

    const upload = async () => {
        if (!name || !file) return;
        setLoading(true);

        const reader = new FileReader();
        reader.onload = async () => {
            try {
                // Parse JSON content if it's a JSON file
                let content = reader.result;
                let parsed = null;

                if (file.name.endsWith('.json')) {
                    try {
                        parsed = JSON.parse(content as string);
                    } catch (e) {
                        console.error("Invalid JSON");
                        alert("Invalid JSON file");
                        setLoading(false);
                        return;
                    }
                }

                // Use the parsed object specifically, e.g. if it has 'atoms' or 'bonds', 
                // we might want to wrap it in 'json_graph' key if the helper expects exact DB schema.
                // However, assuming the backend helper (LibraryAPI.upload) just passes payload
                // and backend expects MoleculeCreate(..., json_graph=..., ...).
                // If the uploaded JSON IS the graph, we wrap it.
                // If the uploaded JSON IS the molecule object (with name etc), we spread it.
                // Let's assume it is the graph export.

                await onUpload({
                    name,
                    json_graph: parsed,
                    previewImage: null
                });

                setName("");
                setFile(null);
            } catch (e) {
                console.error(e);
                alert("Upload failed");
            } finally {
                setLoading(false);
            }
        };

        reader.readAsText(file);
    };

    return (
        <div className="border border-gray-200 p-4 rounded-lg bg-gray-50">
            <div className="font-medium mb-3 flex items-center gap-2">
                <Upload size={16} />
                Upload Molecule
            </div>

            <input
                type="text"
                placeholder="Name"
                className="w-full px-3 py-2 rounded border border-gray-300 mb-3 text-sm"
                value={name}
                onChange={e => setName(e.target.value)}
            />

            <div className="relative mb-3">
                <input
                    type="file"
                    accept=".json"
                    onChange={e => setFile(e.target.files?.[0] || null)}
                    className="block w-full text-sm text-gray-500
              file:mr-4 file:py-2 file:px-4
              file:rounded-full file:border-0
              file:text-sm file:font-semibold
              file:bg-blue-50 file:text-blue-700
              hover:file:bg-blue-100
            "
                />
                {file && <div className="text-xs text-gray-500 mt-1">{file.name}</div>}
            </div>

            <button
                className="w-full py-2 bg-blue-600 text-white rounded font-medium text-sm hover:bg-blue-700 disabled:opacity-50"
                onClick={upload}
                disabled={loading || !file || !name}
            >
                {loading ? "Uploading..." : "Add to Library"}
            </button>
        </div>
    );
}
