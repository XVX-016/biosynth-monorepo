import { apiClient } from '../../api/api';
import type { StudioMessage, StudioMode } from '../../types/studio';

const SKILL_MAP: Record<string, string> = {
    'alina': 'molecule-builder',
    'ionescu': 'property-analyst',
    'chen': 'reaction-planner',
    'solis': 'reaction-simulator',
};

export interface ChatRequest {
    prompt: string;
    mentorId: string;
    mode: StudioMode;
    history?: StudioMessage[];
}

export interface ChatResponse {
    response: string;
    mentorId: string;
    blocked: boolean;
    error?: string;
}

/**
 * Studio AI Gateway
 * 
 * Centralizes all AI interactions for MolForge Studio.
 * Handles mentor-specific skill mapping and error normalization.
 */
export const StudioGateway = {
    chat: async (request: ChatRequest): Promise<ChatResponse> => {
        const { prompt, mentorId } = request;
        const mentorSkill = SKILL_MAP[mentorId] || 'molecule-builder';

        try {
            // In a real production app, we would inject the system prompt here or in the backend
            // For now, we use the existing backend mentor chat endpoint
            const response = await apiClient.post('/api/mentor/chat', {
                prompt: prompt,
                mentor_skill: mentorSkill,
                // We could extend the backend to accept 'mode' and 'system_prompt'
            });

            return {
                response: response.data.response,
                mentorId: mentorId,
                blocked: response.data.blocked || false,
            };
        } catch (error: any) {
            console.error('StudioGateway Chat Error:', error);

            const errorMessage = error.response?.data?.detail || error.message || 'Failed to communicate with mentor';

            return {
                response: '',
                mentorId: mentorId,
                blocked: false,
                error: errorMessage,
            };
        }
    }
};
