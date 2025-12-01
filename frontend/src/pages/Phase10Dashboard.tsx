/**
 * Phase 10 Dashboard
 * 
 * Main dashboard for RL + Generative molecule design.
 */

import React, { useState } from 'react';
import RLWorkflowPanel from '../components/RLWorkflowPanel';
import TopCandidatesPanel from '../components/TopCandidatesPanel';
import RewardVisualization from '../components/RewardVisualization';
import type { RunLoopResponse, TopCandidate } from '../api/phase10';

export default function Phase10Dashboard() {
  const [refreshTrigger, setRefreshTrigger] = useState(0);
  const [selectedCandidate, setSelectedCandidate] = useState<TopCandidate | null>(null);
  const [workflowResults, setWorkflowResults] = useState<RunLoopResponse | null>(null);

  const handleWorkflowComplete = (results: RunLoopResponse) => {
    setWorkflowResults(results);
    setRefreshTrigger(prev => prev + 1);
  };

  const handleSelectCandidate = (candidate: TopCandidate) => {
    setSelectedCandidate(candidate);
    // TODO: Open molecule in 3D viewer or lab
    console.log('Selected candidate:', candidate);
  };

  return (
    <div className="min-h-screen bg-gray-50 p-6">
      <div className="max-w-7xl mx-auto">
        <div className="mb-6">
          <h1 className="text-3xl font-bold text-gray-900">Phase 10: RL + Generative Design</h1>
          <p className="text-gray-600 mt-2">
            Reinforcement learning and generative model-based molecule optimization
          </p>
        </div>

        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          {/* Left Column: Workflow Control */}
          <div className="lg:col-span-1">
            <RLWorkflowPanel onWorkflowComplete={handleWorkflowComplete} />
          </div>

          {/* Middle Column: Top Candidates */}
          <div className="lg:col-span-1">
            <TopCandidatesPanel
              refreshTrigger={refreshTrigger}
              onSelectCandidate={handleSelectCandidate}
            />
          </div>

          {/* Right Column: Visualization */}
          <div className="lg:col-span-1">
            <RewardVisualization refreshTrigger={refreshTrigger} />
          </div>
        </div>

        {/* Workflow Results Summary */}
        {workflowResults && (
          <div className="mt-6 bg-white rounded-lg shadow p-6">
            <h2 className="text-xl font-bold mb-4">Latest Workflow Results</h2>
            <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
              <div>
                <div className="text-sm text-gray-500">Iterations</div>
                <div className="text-2xl font-bold">
                  {workflowResults.iterations_completed}
                </div>
              </div>
              <div>
                <div className="text-sm text-gray-500">Top Candidates</div>
                <div className="text-2xl font-bold">
                  {workflowResults.top_candidates.length}
                </div>
              </div>
              <div>
                <div className="text-sm text-gray-500">Best Reward</div>
                <div className="text-2xl font-bold text-green-600">
                  {workflowResults.top_candidates[0]?.reward.toFixed(4) || '0.0000'}
                </div>
              </div>
              <div>
                <div className="text-sm text-gray-500">Total Molecules</div>
                <div className="text-2xl font-bold">
                  {workflowResults.statistics?.total || 0}
                </div>
              </div>
            </div>
          </div>
        )}

        {/* Selected Candidate Details */}
        {selectedCandidate && (
          <div className="mt-6 bg-white rounded-lg shadow p-6">
            <h2 className="text-xl font-bold mb-4">Selected Candidate</h2>
            <div className="grid grid-cols-2 gap-4">
              <div>
                <div className="text-sm text-gray-500">SMILES</div>
                <div className="font-mono text-lg">{selectedCandidate.smiles}</div>
              </div>
              <div>
                <div className="text-sm text-gray-500">Reward</div>
                <div className="text-2xl font-bold text-green-600">
                  {selectedCandidate.reward.toFixed(4)}
                </div>
              </div>
              <div>
                <div className="text-sm text-gray-500">Iteration</div>
                <div className="text-lg">{selectedCandidate.iteration}</div>
              </div>
              <div>
                <div className="text-sm text-gray-500">Method</div>
                <div className="text-lg">{selectedCandidate.method}</div>
              </div>
            </div>
            {Object.keys(selectedCandidate.properties).length > 0 && (
              <div className="mt-4">
                <div className="text-sm text-gray-500 mb-2">Properties</div>
                <div className="grid grid-cols-3 gap-2">
                  {Object.entries(selectedCandidate.properties).map(([key, value]) => (
                    <div key={key} className="bg-gray-50 p-2 rounded">
                      <div className="text-xs text-gray-500">{key}</div>
                      <div className="font-medium">
                        {typeof value === 'number' ? value.toFixed(3) : value}
                      </div>
                    </div>
                  ))}
                </div>
              </div>
            )}
          </div>
        )}
      </div>
    </div>
  );
}

