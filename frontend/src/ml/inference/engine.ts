import type { Molecule } from "../../chemcore/graph/molecule";
import { extractFeatures } from "../features/extractor";
import type { MLModel, PredictionResult } from "../predictors/types";

export class InferenceEngine {
    private models: Map<string, MLModel> = new Map();

    register(model: MLModel) {
        this.models.set(model.id, model);
    }

    async run(molecule: Molecule, modelId: string): Promise<PredictionResult> {
        const model = this.models.get(modelId);
        if (!model) throw new Error(`Model ${modelId} not found`);

        const features = extractFeatures(molecule);
        return model.predict(features);
    }

    listModels() {
        return Array.from(this.models.values());
    }
}

export const inferenceEngine = new InferenceEngine();
