import { apiClient } from './api';
import type { MLRequest, OptimizeResponse, BondOrderResponse, AutoBondResponse } from "../types/ml";

export const MLAPI = {
    optimize: async (payload: MLRequest): Promise<OptimizeResponse> => {
        const r = await apiClient.post("/ml/optimize", payload);
        return r.data;
    },
    bondOrder: async (payload: MLRequest): Promise<BondOrderResponse> => {
        const r = await apiClient.post("/ml/bond-order", payload);
        return r.data;
    },
    autobond: async (payload: MLRequest): Promise<AutoBondResponse> => {
        const r = await apiClient.post("/ml/autobond", payload);
        return r.data;
    },
    predictProperties: async (payload: MLRequest) => {
        const r = await apiClient.post("/lab/run-ml-prediction", payload);
        return r.data;
    }
};
