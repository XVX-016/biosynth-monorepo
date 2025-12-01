/**
 * Top Candidates Panel
 * 
 * Displays top candidates by reward score.
 */

import React, { useEffect, useState } from 'react';
import { getTopCandidates, type TopCandidate } from '../api/phase10';

interface TopCandidatesPanelProps {
  refreshTrigger?: number;
  onSelectCandidate?: (candidate: TopCandidate) => void;
}

export default function TopCandidatesPanel({
  refreshTrigger,
  onSelectCandidate,
}: TopCandidatesPanelProps) {
  const [candidates, setCandidates] = useState<TopCandidate[]>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const loadCandidates = async () => {
    setLoading(true);
    setError(null);

    try {
      const top = await getTopCandidates(10);
      setCandidates(top);
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to load candidates');
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    loadCandidates();
  }, [refreshTrigger]);

  if (loading) {
    return (
      <div className="bg-white rounded-lg shadow p-6">
        <h2 className="text-xl font-bold mb-4">Top Candidates</h2>
        <div className="text-center py-8">Loading...</div>
      </div>
    );
  }

  if (error) {
    return (
      <div className="bg-white rounded-lg shadow p-6">
        <h2 className="text-xl font-bold mb-4">Top Candidates</h2>
        <div className="bg-red-50 text-red-700 px-4 py-2 rounded">{error}</div>
      </div>
    );
  }

  return (
    <div className="bg-white rounded-lg shadow p-6">
      <div className="flex justify-between items-center mb-4">
        <h2 className="text-xl font-bold">Top Candidates</h2>
        <button
          onClick={loadCandidates}
          className="text-sm text-blue-600 hover:text-blue-700"
        >
          Refresh
        </button>
      </div>

      {candidates.length === 0 ? (
        <div className="text-center py-8 text-gray-500">
          No candidates yet. Run a workflow loop to generate molecules.
        </div>
      ) : (
        <div className="space-y-3">
          {candidates.map((candidate, idx) => (
            <div
              key={idx}
              className="border rounded p-3 hover:bg-gray-50 cursor-pointer transition"
              onClick={() => onSelectCandidate?.(candidate)}
            >
              <div className="flex justify-between items-start mb-2">
                <div className="flex-1">
                  <div className="font-mono text-sm text-gray-700 mb-1">
                    {candidate.smiles}
                  </div>
                  <div className="text-xs text-gray-500">
                    Iteration {candidate.iteration} • {candidate.method}
                  </div>
                </div>
                <div className="text-right">
                  <div className="text-lg font-bold text-green-600">
                    {candidate.reward.toFixed(4)}
                  </div>
                  <div className="text-xs text-gray-500">Reward</div>
                </div>
              </div>

              {Object.keys(candidate.properties).length > 0 && (
                <div className="mt-2 pt-2 border-t text-xs">
                  <div className="flex flex-wrap gap-2">
                    {Object.entries(candidate.properties).map(([key, value]) => (
                      <span key={key} className="text-gray-600">
                        {key}: {typeof value === 'number' ? value.toFixed(2) : value}
                      </span>
                    ))}
                  </div>
                </div>
              )}
            </div>
          ))}
        </div>
      )}
    </div>
  );
}

