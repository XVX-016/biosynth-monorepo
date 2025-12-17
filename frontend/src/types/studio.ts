/**
 * MolForge Studio - Architectural Contracts
 */

export type StudioMode = 'design' | 'optimize' | 'simulate';

export interface StudioMentor {
    id: string;
    name: string;
    role: string;
    personality: string;
    allowedDomains: string[];
    systemPrompt: string;
    avatarModelPath: string;
    accentColor: string;
    backgroundGradient: string;
}

export interface StudioMessage {
    id: string;
    role: 'user' | 'mentor';
    content: string;
    timestamp: Date;
    metadata?: {
        moleculeId?: string;
        smiles?: string;
        mode?: StudioMode;
    };
}

export interface StudioSession {
    id: string;
    activeMentorId: string;
    currentMode: StudioMode;
    messages: StudioMessage[];
    lastContext?: {
        molfile?: string;
        smiles?: string;
    };
}
