import { apiClient } from "./api"; // Use centralized client
import type { Molecule } from "../types/molecule";

export const LibraryAPI = {
    list: async (params?: { page?: number; limit?: number }) => {
        const p = params?.page || 1;
        const l = params?.limit || 9; // Default to 9 for 3x3 grid
        const offset = (p - 1) * l;

        const res = await apiClient.get(`/molecules/list?limit=${l}&offset=${offset}`);

        // Adapt backend response to frontend Molecule interface
        return res.data.map((m: any) => ({
            id: String(m.id),
            name: m.name,
            formula: m.formula || "",
            // Assume backend sends raw base64, prepend data URI scheme if needed. 
            // If backend sends null, empty string.
            previewImage: m.thumbnail_b64 ? (m.thumbnail_b64.startsWith('data:') ? m.thumbnail_b64 : `data:image/png;base64,${m.thumbnail_b64}`) : "",
            createdAt: m.created_at,
            updatedAt: m.created_at,
            atoms: [], // List view doesn't need full structure
            bonds: [],
            isValid: true,
            qualityScore: 100
        }));
    },

    upload: async (payload: any) => {
        const res = await apiClient.post("/molecules/save", payload);
        return res.data;
    },

    delete: async (id: string) => {
        const res = await apiClient.delete(`/molecules/${id}`);
        return res.data;
    },

    get: async (id: string) => {
        const res = await apiClient.get(`/molecules/${id}`);
        return res.data;
    },

    export: async (id: string) => {
        const res = await apiClient.post(`/molecules/${id}/export`);
        return res.data;
    },
};
