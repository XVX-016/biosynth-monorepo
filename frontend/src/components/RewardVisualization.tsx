/**
 * Reward Visualization
 * 
 * Charts and graphs for reward progression.
 */

import React, { useEffect, useState } from 'react';
import { getIterationLogs, getStatistics } from '../api/phase10';

interface RewardVisualizationProps {
  refreshTrigger?: number;
}

export default function RewardVisualization({ refreshTrigger }: RewardVisualizationProps) {
  const [logs, setLogs] = useState<Array<Record<string, any>>>([]);
  const [stats, setStats] = useState<Record<string, any>>({});
  const [loading, setLoading] = useState(false);

  const loadData = async () => {
    setLoading(true);
    try {
      const [logsData, statsData] = await Promise.all([
        getIterationLogs(),
        getStatistics(),
      ]);
      setLogs(logsData);
      setStats(statsData);
    } catch (err) {
      console.error('Failed to load visualization data:', err);
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    loadData();
  }, [refreshTrigger]);

  if (loading) {
    return (
      <div className="bg-white rounded-lg shadow p-6">
        <h2 className="text-xl font-bold mb-4">Reward Progression</h2>
        <div className="text-center py-8">Loading...</div>
      </div>
    );
  }

  if (logs.length === 0) {
    return (
      <div className="bg-white rounded-lg shadow p-6">
        <h2 className="text-xl font-bold mb-4">Reward Progression</h2>
        <div className="text-center py-8 text-gray-500">
          No iteration data yet. Run a workflow loop to see progress.
        </div>
      </div>
    );
  }

  // Extract reward data for chart
  const avgRewards = logs.map(log => log.avg_reward || 0);
  const maxRewards = logs.map(log => log.max_reward || 0);
  const iterations = logs.map(log => log.iteration || 0);

  // Simple bar chart visualization
  const maxReward = Math.max(...maxRewards, 1);

  return (
    <div className="bg-white rounded-lg shadow p-6">
      <h2 className="text-xl font-bold mb-4">Reward Progression</h2>

      {/* Statistics */}
      {stats && (
        <div className="grid grid-cols-3 gap-4 mb-6">
          <div className="text-center">
            <div className="text-2xl font-bold text-blue-600">
              {stats.total || 0}
            </div>
            <div className="text-sm text-gray-500">Total Molecules</div>
          </div>
          <div className="text-center">
            <div className="text-2xl font-bold text-green-600">
              {stats.max_reward?.toFixed(4) || '0.0000'}
            </div>
            <div className="text-sm text-gray-500">Max Reward</div>
          </div>
          <div className="text-center">
            <div className="text-2xl font-bold text-purple-600">
              {stats.avg_reward?.toFixed(4) || '0.0000'}
            </div>
            <div className="text-sm text-gray-500">Avg Reward</div>
          </div>
        </div>
      )}

      {/* Simple bar chart */}
      <div className="space-y-2">
        <div className="text-sm font-medium mb-2">Reward by Iteration</div>
        <div className="space-y-1">
          {iterations.map((iter, idx) => {
            const avgWidth = (avgRewards[idx] / maxReward) * 100;
            const maxWidth = (maxRewards[idx] / maxReward) * 100;

            return (
              <div key={iter} className="flex items-center gap-2">
                <div className="w-12 text-xs text-gray-600">#{iter}</div>
                <div className="flex-1 relative h-6 bg-gray-100 rounded">
                  <div
                    className="absolute h-full bg-blue-400 rounded"
                    style={{ width: `${avgWidth}%` }}
                    title={`Avg: ${avgRewards[idx].toFixed(4)}`}
                  />
                  <div
                    className="absolute h-full bg-green-500 rounded opacity-75"
                    style={{ width: `${maxWidth}%` }}
                    title={`Max: ${maxRewards[idx].toFixed(4)}`}
                  />
                </div>
                <div className="w-20 text-xs text-gray-600 text-right">
                  {maxRewards[idx].toFixed(3)}
                </div>
              </div>
            );
          })}
        </div>
      </div>

      {/* Legend */}
      <div className="flex gap-4 mt-4 text-xs">
        <div className="flex items-center gap-2">
          <div className="w-4 h-4 bg-blue-400 rounded" />
          <span>Average Reward</span>
        </div>
        <div className="flex items-center gap-2">
          <div className="w-4 h-4 bg-green-500 rounded" />
          <span>Max Reward</span>
        </div>
      </div>
    </div>
  );
}

