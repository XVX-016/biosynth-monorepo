/**
 * RL Workflow Panel
 * 
 * Controls for running RL workflow loops.
 */

import React, { useState } from 'react';
import { runWorkflowLoop, type RunLoopRequest } from '../api/phase10';

interface RLWorkflowPanelProps {
  onWorkflowComplete?: (results: any) => void;
}

export default function RLWorkflowPanel({ onWorkflowComplete }: RLWorkflowPanelProps) {
  const [maxIterations, setMaxIterations] = useState(10);
  const [batchSize, setBatchSize] = useState(32);
  const [useGenerative, setUseGenerative] = useState(false);
  const [seedSmiles, setSeedSmiles] = useState('');
  const [running, setRunning] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const handleRun = async () => {
    setRunning(true);
    setError(null);

    try {
      const request: RunLoopRequest = {
        max_iterations: maxIterations,
        batch_size: batchSize,
        use_generative: useGenerative,
        seed_smiles: seedSmiles
          ? seedSmiles.split(',').map(s => s.trim()).filter(Boolean)
          : undefined,
      };

      const results = await runWorkflowLoop(request);
      onWorkflowComplete?.(results);
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Unknown error');
    } finally {
      setRunning(false);
    }
  };

  return (
    <div className="bg-white rounded-lg shadow p-6">
      <h2 className="text-xl font-bold mb-4">RL Workflow Control</h2>

      <div className="space-y-4">
        <div>
          <label className="block text-sm font-medium mb-1">
            Max Iterations
          </label>
          <input
            type="number"
            value={maxIterations}
            onChange={(e) => setMaxIterations(parseInt(e.target.value) || 10)}
            className="w-full px-3 py-2 border rounded"
            min="1"
            max="100"
            disabled={running}
          />
        </div>

        <div>
          <label className="block text-sm font-medium mb-1">
            Batch Size
          </label>
          <input
            type="number"
            value={batchSize}
            onChange={(e) => setBatchSize(parseInt(e.target.value) || 32)}
            className="w-full px-3 py-2 border rounded"
            min="1"
            max="1000"
            disabled={running}
          />
        </div>

        <div>
          <label className="flex items-center gap-2">
            <input
              type="checkbox"
              checked={useGenerative}
              onChange={(e) => setUseGenerative(e.target.checked)}
              disabled={running}
            />
            <span className="text-sm font-medium">Use Generative Agent</span>
          </label>
        </div>

        <div>
          <label className="block text-sm font-medium mb-1">
            Seed SMILES (comma-separated, optional)
          </label>
          <input
            type="text"
            value={seedSmiles}
            onChange={(e) => setSeedSmiles(e.target.value)}
            className="w-full px-3 py-2 border rounded"
            placeholder="CCO, CCCO"
            disabled={running}
          />
        </div>

        {error && (
          <div className="bg-red-50 text-red-700 px-4 py-2 rounded">
            {error}
          </div>
        )}

        <button
          onClick={handleRun}
          disabled={running}
          className="w-full bg-blue-600 text-white px-4 py-2 rounded hover:bg-blue-700 disabled:opacity-50 disabled:cursor-not-allowed"
        >
          {running ? 'Running...' : 'Run Workflow Loop'}
        </button>
      </div>
    </div>
  );
}

