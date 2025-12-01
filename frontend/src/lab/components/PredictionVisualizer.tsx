/**
 * PredictionVisualizer - Visualize predictions with heatmaps and contribution maps
 */

import React, { useMemo } from 'react';

export interface PredictionVisualizerProps {
  predictions: Record<string, number>;
  contributions?: Record<string, number>; // Atom/bond contributions
  attentions?: {
    edgeAttentions: number[]; // Per-edge attention weights
    edgeIndex: number[][]; // Edge indices [2, E]
    nodeImportance?: number[]; // Per-node importance
  };
  showHeatmap?: boolean;
  showContributions?: boolean;
  showAttention?: boolean;
}

export default function PredictionVisualizer({
  predictions,
  contributions,
  attentions,
  showHeatmap = true,
  showContributions = false,
  showAttention = false,
}: PredictionVisualizerProps) {
  // Normalize predictions for visualization
  const normalized = useMemo(() => {
    const values = Object.values(predictions);
    const min = Math.min(...values);
    const max = Math.max(...values);
    const range = max - min || 1;

    return Object.entries(predictions).reduce((acc, [key, value]) => {
      acc[key] = {
        value,
        normalized: (value - min) / range,
      };
      return acc;
    }, {} as Record<string, { value: number; normalized: number }>);
  }, [predictions]);

  // Color scale (blue to red)
  const getColor = (normalized: number): string => {
    if (normalized < 0.5) {
      // Blue to white
      const intensity = normalized * 2;
      const r = Math.floor(255 * intensity);
      const g = Math.floor(255 * intensity);
      const b = 255;
      return `rgb(${r}, ${g}, ${b})`;
    } else {
      // White to red
      const intensity = (normalized - 0.5) * 2;
      const r = 255;
      const g = Math.floor(255 * (1 - intensity));
      const b = Math.floor(255 * (1 - intensity));
      return `rgb(${r}, ${g}, ${b})`;
    }
  };

  // Format value
  const formatValue = (value: number): string => {
    if (Math.abs(value) < 0.01) {
      return value.toExponential(2);
    }
    return value.toFixed(3);
  };

  return (
    <div className="space-y-3">
      {/* Predictions List */}
      <div className="space-y-2">
        <div className="text-xs font-semibold opacity-70">Predictions</div>
        {Object.entries(normalized).map(([property, { value, normalized: norm }]) => (
          <div key={property} className="bg-neutral-700 rounded p-2">
            <div className="flex items-center justify-between mb-1">
              <span className="text-xs font-medium text-neutral-200">{property}</span>
              <span className="text-xs font-mono text-neutral-300">{formatValue(value)}</span>
            </div>
            {showHeatmap && (
              <div className="mt-1">
                <div
                  className="h-2 rounded"
                  style={{
                    width: `${norm * 100}%`,
                    backgroundColor: getColor(norm),
                  }}
                />
              </div>
            )}
          </div>
        ))}
      </div>

      {/* Heatmap Visualization */}
      {showHeatmap && (
        <div className="mt-3">
          <div className="text-xs font-semibold opacity-70 mb-2">Heatmap</div>
          <div className="grid grid-cols-2 gap-2">
            {Object.entries(normalized).map(([property, { normalized: norm }]) => (
              <div
                key={property}
                className="p-2 rounded text-center"
                style={{
                  backgroundColor: getColor(norm),
                  color: norm > 0.5 ? '#fff' : '#000',
                }}
              >
                <div className="text-xs font-medium">{property}</div>
                <div className="text-xs opacity-80">{norm.toFixed(2)}</div>
              </div>
            ))}
          </div>
        </div>
      )}

      {/* Contribution Map */}
      {showContributions && contributions && (
        <div className="mt-3">
          <div className="text-xs font-semibold opacity-70 mb-2">Contributions</div>
          <div className="space-y-1">
            {Object.entries(contributions)
              .sort(([, a], [, b]) => Math.abs(b) - Math.abs(a))
              .slice(0, 10)
              .map(([atomId, contribution]) => (
                <div key={atomId} className="flex items-center justify-between text-xs">
                  <span className="opacity-60">{atomId.slice(0, 8)}</span>
                  <span
                    className={`font-mono ${
                      contribution > 0 ? 'text-green-400' : 'text-red-400'
                    }`}
                  >
                    {contribution > 0 ? '+' : ''}
                    {formatValue(contribution)}
                  </span>
                </div>
              ))}
          </div>
        </div>
      )}

      {/* Attention Visualization */}
      {showAttention && attentions && (
        <div className="mt-3">
          <div className="text-xs font-semibold opacity-70 mb-2">Attention Weights</div>
          <div className="text-xs opacity-60 mb-2">
            Showing attention for {attentions.edgeAttentions.length} bonds
          </div>
          {attentions.nodeImportance && (
            <div className="mb-2">
              <div className="text-xs opacity-60 mb-1">Node Importance (Top 10)</div>
              <div className="space-y-1">
                {attentions.nodeImportance
                  .map((importance, idx) => ({ idx, importance }))
                  .sort((a, b) => b.importance - a.importance)
                  .slice(0, 10)
                  .map(({ idx, importance }) => (
                    <div key={idx} className="flex items-center justify-between text-xs">
                      <span className="opacity-60">Atom {idx}</span>
                      <div className="flex items-center gap-2">
                        <div
                          className="h-2 w-16 rounded"
                          style={{
                            backgroundColor: getColor(importance),
                          }}
                        />
                        <span className="font-mono">{importance.toFixed(3)}</span>
                      </div>
                    </div>
                  ))}
              </div>
            </div>
          )}
          <div className="text-xs opacity-60 mt-2">
            Attention weights indicate which bonds/atoms the model focuses on for predictions.
            Higher values = more important.
          </div>
        </div>
      )}

      {/* Color Legend */}
      {showHeatmap && (
        <div className="mt-3 pt-3 border-t border-neutral-700">
          <div className="text-xs opacity-60 mb-1">Color Scale</div>
          <div className="flex items-center gap-1">
            <div className="text-xs">Low</div>
            <div
              className="flex-1 h-3 rounded"
              style={{
                background: 'linear-gradient(to right, rgb(0,0,255), rgb(255,255,255), rgb(255,0,0))',
              }}
            />
            <div className="text-xs">High</div>
          </div>
        </div>
      )}
    </div>
  );
}

