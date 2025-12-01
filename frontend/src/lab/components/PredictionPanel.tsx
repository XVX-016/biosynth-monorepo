/**
 * PredictionPanel - Interactive prediction UI panel
 */

import React, { useState, useEffect } from 'react';
import { useLab } from '../hooks/useLab';
import { predict } from '../api/predict';
import { predictWithAttention } from '../api/predictWithAttention';
import { modelRegistry } from '../ml/ModelRegistry';
import type { ModelInfo } from '../ml/ModelRegistry';
import PredictionVisualizer from './PredictionVisualizer';
import { useAttentionWeights } from './AttentionVisualizer';

export default function PredictionPanel() {
  const { currentMolecule } = useLab();
  const [models, setModels] = useState<ModelInfo[]>([]);
  const [selectedModel, setSelectedModel] = useState<string>('');
  const [predicting, setPredicting] = useState(false);
  const [predictions, setPredictions] = useState<Record<string, number> | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [properties, setProperties] = useState<string[]>(['logP', 'toxicity', 'solubility']);
  const [attentions, setAttentions] = useState<{
    edgeAttentions: number[];
    edgeIndex: number[][];
    nodeImportance?: number[];
  } | null>(null);
  const [useAttention, setUseAttention] = useState(false);
  const { setAttentionWeights } = useLab();

  useEffect(() => {
    // Load models from registry
    modelRegistry.initDefaults().then(() => {
      const modelList = modelRegistry.listModels();
      setModels(modelList);
      if (modelList.length > 0 && !selectedModel) {
        const defaultModel = modelRegistry.getDefault();
        setSelectedModel(defaultModel?.id || modelList[0].id);
      }
    });
  }, []);

  const handlePredict = async () => {
    if (!currentMolecule || currentMolecule.atoms.size === 0) {
      setError('No molecule loaded');
      return;
    }

    if (!selectedModel) {
      setError('No model selected');
      return;
    }

    setPredicting(true);
    setError(null);
    setPredictions(null);
    setAttentions(null);

    try {
      if (useAttention && selectedModelInfo?.type === 'gnn') {
        // Use attention-based prediction
        const result = await predictWithAttention({
          molecule: currentMolecule,
          modelId: selectedModel,
          returnAttention: true,
        });

        // Convert single prediction to property dict (simplified)
        setPredictions({
          prediction: Array.isArray(result.prediction)
            ? result.prediction[0]
            : result.prediction,
        });

        // Normalize attentions
        const normalized = result.attentions[result.attentions.length - 1] || [];
        const min = Math.min(...normalized);
        const max = Math.max(...normalized);
        const range = max - min || 1;
        const normalizedAttentions = normalized.map((att) => (att - min) / range);

        const attentionData = {
          edgeAttentions: normalizedAttentions,
          edgeIndex: result.edgeIndex,
          nodeImportance: result.nodeImportance,
        };
        
        setAttentions(attentionData);
        
        // Map attention weights to bond/atom IDs and store in LabStore
        if (currentMolecule && attentionData.edgeIndex && attentionData.edgeAttentions) {
          const { mapAttentionWeights } = await import('../utils/mapAttention');
          
          const atoms = Array.from(currentMolecule.atoms.values());
          const bonds = Array.from(currentMolecule.bonds.values());
          
          const { bondAttentions, atomImportance } = mapAttentionWeights(
            attentionData.edgeIndex,
            attentionData.edgeAttentions,
            attentionData.nodeImportance,
            atoms,
            bonds
          );
          
          setAttentionWeights({ bondAttentions, atomImportance });
        }
      } else {
        // Use standard prediction
        const result = await predict({
          molecule: currentMolecule,
          modelId: selectedModel,
          properties,
        });

        setPredictions(result.predictions);
      }
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Prediction failed');
    } finally {
      setPredicting(false);
    }
  };

  const selectedModelInfo = models.find((m) => m.id === selectedModel);

  return (
    <div className="p-3 border-t border-neutral-700">
      <h2 className="text-sm font-semibold mb-2 opacity-70">Predictions</h2>

      {/* Model Selection */}
      <div className="mb-3">
        <label className="text-xs opacity-60 mb-1 block">Model</label>
        <select
          value={selectedModel}
          onChange={(e) => setSelectedModel(e.target.value)}
          className="w-full px-2 py-1 rounded bg-neutral-700 text-neutral-200 text-xs"
        >
          {models.map((model) => (
            <option key={model.id} value={model.id}>
              {model.name} ({model.type})
            </option>
          ))}
        </select>
        {selectedModelInfo && (
          <div className="text-xs opacity-60 mt-1">
            {selectedModelInfo.description}
          </div>
        )}
      </div>

      {/* Attention Toggle */}
      {selectedModelInfo?.type === 'gnn' && (
        <div className="mb-3">
          <label className="flex items-center text-xs">
            <input
              type="checkbox"
              checked={useAttention}
              onChange={(e) => setUseAttention(e.target.checked)}
              className="mr-2"
            />
            Show attention weights (explainability)
          </label>
        </div>
      )}

      {/* Properties Selection */}
      <div className="mb-3">
        <label className="text-xs opacity-60 mb-1 block">Properties</label>
        <div className="space-y-1">
          {['logP', 'toxicity', 'solubility', 'molecularWeight', 'polarSurfaceArea'].map(
            (prop) => (
              <label key={prop} className="flex items-center text-xs">
                <input
                  type="checkbox"
                  checked={properties.includes(prop)}
                  onChange={(e) => {
                    if (e.target.checked) {
                      setProperties([...properties, prop]);
                    } else {
                      setProperties(properties.filter((p) => p !== prop));
                    }
                  }}
                  className="mr-2"
                />
                {prop}
              </label>
            )
          )}
        </div>
      </div>

      {/* Predict Button */}
      <button
        onClick={handlePredict}
        disabled={predicting || !currentMolecule || currentMolecule.atoms.size === 0}
        className="w-full px-3 py-2 rounded bg-blue-600 text-white text-sm font-medium hover:bg-blue-700 disabled:opacity-50 disabled:cursor-not-allowed mb-3"
      >
        {predicting ? 'Predicting...' : 'Predict Properties'}
      </button>

      {/* Error Display */}
      {error && (
        <div className="mb-3 p-2 rounded bg-red-900/50 text-red-300 text-xs">
          {error}
        </div>
      )}

      {/* Predictions Display */}
      {predictions && (
        <div className="mt-3">
          <PredictionVisualizer
            predictions={predictions}
            attentions={attentions || undefined}
            showAttention={useAttention && !!attentions}
          />
        </div>
      )}

      {/* No Molecule Warning */}
      {(!currentMolecule || currentMolecule.atoms.size === 0) && (
        <div className="text-xs opacity-60 text-center py-2">
          Load a molecule to make predictions
        </div>
      )}
    </div>
  );
}

