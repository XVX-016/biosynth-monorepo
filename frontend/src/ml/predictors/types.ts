import type { MoleculeFeatures } from "../features/extractor";

export interface PredictionResult {
    value: number | boolean;
    confidence?: number;
    units?: string;
}

export interface MLModel {
    id: string;
    name: string;
    description: string;
    type: 'regression' | 'classification';
    version: string;
    predict(features: MoleculeFeatures): Promise<PredictionResult>;
}
