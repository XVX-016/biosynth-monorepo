// src/api/ml.ts
import axios from "axios";

// Create axios instance if not existing, or reuse
const api = axios.create({
    baseURL: "/api" // Adjust if needed
});

import type { MLRequest, OptimizeResponse, BondOrderResponse, AutoBondResponse } from "../types/ml";

export const MLAPI = {
    optimize: async (payload: MLRequest): Promise<OptimizeResponse> => {
        // Note: app.py mounts ml router at /api/ml, and if we removed prefix in router, it is /api/ml/optimize
        const r = await api.post("/ml/optimize", payload);
        return r.data;
    },
    bondOrder: async (payload: MLRequest): Promise<BondOrderResponse> => {
        const r = await api.post("/ml/bond-order", payload);
        return r.data;
    },
    autobond: async (payload: MLRequest): Promise<AutoBondResponse> => {
        const r = await api.post("/ml/autobond", payload);
        return r.data;
    },
    predictProperties: async (payload: MLRequest) => {
        const r = await api.post("/ml/predict-properties", payload);
        return r.data;
    }
};
