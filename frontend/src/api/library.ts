import axios from "axios";

// Helper to get base URL if needed, or rely on proxy/relative path
const API_URL = "/api";

export const LibraryAPI = {
    list: async () => {
        const res = await axios.get("/molecules/list");
        return res.data;
    },

    upload: async (payload: any) => {
        const res = await axios.post("/molecules/save", payload);
        return res.data;
    },

    delete: async (id: string) => {
        const res = await axios.delete(`/molecules/${id}`);
        return res.data;
    },

    get: async (id: string) => {
        const res = await axios.get(`/molecules/${id}`);
        return res.data;
    },

    export: async (id: string) => {
        const res = await axios.post(`/molecules/${id}/export`);
        return res.data;
    },
};
