/**
 * Phase 10 API Client
 * 
 * Client for RL + Generative molecule design endpoints.
 */

const API_BASE = import.meta.env.VITE_API_URL || 'http://localhost:8000';

export interface GenerateRequest {
  n: number;
  method: 'rl' | 'generative';
  seed_smiles?: string[];
}

export interface GenerateResponse {
  molecules: string[];
  method: string;
  count: number;
}

export interface EvaluateRequest {
  smiles: string;
  compute_ml?: boolean;
  compute_screening?: boolean;
  compute_qm?: boolean;
  compute_md?: boolean;
}

export interface EvaluateResponse {
  smiles: string;
  reward: number;
  ml_predictions: Record<string, number>;
  screening_results: Record<string, any>;
  qm_results: Record<string, number>;
  md_results: Record<string, number>;
}

export interface RunLoopRequest {
  max_iterations: number;
  batch_size?: number;
  use_generative?: boolean;
  seed_smiles?: string[];
}

export interface RunLoopResponse {
  iterations_completed: number;
  top_candidates: Array<{
    smiles: string;
    reward: number;
    properties: Record<string, number>;
    iteration: number;
  }>;
  statistics: Record<string, any>;
  iteration_logs: Array<Record<string, any>>;
}

export interface TopCandidate {
  smiles: string;
  reward: number;
  properties: Record<string, number>;
  iteration: number;
  method: string;
}

/**
 * Generate molecules using RL or generative agent.
 */
export async function generateMolecules(
  request: GenerateRequest
): Promise<GenerateResponse> {
  const response = await fetch(`${API_BASE}/api/phase10/generate`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(request),
  });

  if (!response.ok) {
    throw new Error(`Generation failed: ${response.statusText}`);
  }

  return response.json();
}

/**
 * Evaluate a single molecule.
 */
export async function evaluateMolecule(
  request: EvaluateRequest
): Promise<EvaluateResponse> {
  const response = await fetch(`${API_BASE}/api/phase10/evaluate`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(request),
  });

  if (!response.ok) {
    throw new Error(`Evaluation failed: ${response.statusText}`);
  }

  return response.json();
}

/**
 * Run the RL workflow loop.
 */
export async function runWorkflowLoop(
  request: RunLoopRequest
): Promise<RunLoopResponse> {
  const response = await fetch(`${API_BASE}/api/phase10/run_loop`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(request),
  });

  if (!response.ok) {
    throw new Error(`Workflow loop failed: ${response.statusText}`);
  }

  return response.json();
}

/**
 * Get top candidates by reward.
 */
export async function getTopCandidates(n: number = 10): Promise<TopCandidate[]> {
  const response = await fetch(`${API_BASE}/api/phase10/top_candidates?n=${n}`);

  if (!response.ok) {
    throw new Error(`Failed to get top candidates: ${response.statusText}`);
  }

  const data = await response.json();
  return data.candidates;
}

/**
 * Get workflow statistics.
 */
export async function getStatistics(): Promise<Record<string, any>> {
  const response = await fetch(`${API_BASE}/api/phase10/statistics`);

  if (!response.ok) {
    throw new Error(`Failed to get statistics: ${response.statusText}`);
  }

  return response.json();
}

/**
 * Get iteration logs.
 */
export async function getIterationLogs(): Promise<Array<Record<string, any>>> {
  const response = await fetch(`${API_BASE}/api/phase10/iteration_logs`);

  if (!response.ok) {
    throw new Error(`Failed to get iteration logs: ${response.statusText}`);
  }

  const data = await response.json();
  return data.logs;
}

