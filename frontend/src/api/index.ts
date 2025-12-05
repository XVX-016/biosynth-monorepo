const BASE = import.meta.env.VITE_BACKEND_URL || "http://localhost:8000";

export const api = {
    predict: (payload) =>
        fetch(`${BASE}/api/predict/property`, {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify(payload),
        }).then(r => r.json()),

    predictBatch: (payload) =>
        fetch(`${BASE}/api/predict/batch`, {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify(payload),
        }).then(r => r.json()),

    relax: (molfile) =>
        fetch(`${BASE}/api/molecule/relax`, {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ molfile }),
        }).then(r => r.json()),

    generate3D: (smiles) =>
        fetch(`${BASE}/api/molecule/generate-3d`, {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ smiles }),
        }).then(r => r.json()),

    train: () =>
        fetch(`${BASE}/api/models/train`, { method: "POST" }),

    listModels: () =>
        fetch(`${BASE}/api/models/list`).then(r => r.json()),
};
