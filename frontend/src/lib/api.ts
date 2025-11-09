import axios from 'axios'

const API_URL = import.meta.env.VITE_API_URL || 'http://localhost:8000'

const apiClient = axios.create({
  baseURL: API_URL,
  headers: {
    'Content-Type': 'application/json',
  },
})

export interface PredictResponse {
  properties: {
    stability: number
    toxicity: number
    solubility: number
    bioavailability: number
    novelty: number
  }
}

export interface GenerateResponse {
  smiles: string
}

/**
 * Predict molecular properties from SMILES string
 */
export async function predict(smiles: string): Promise<PredictResponse> {
  const response = await apiClient.post<PredictResponse>('/predict', {
    smiles,
  })
  return response.data
}

/**
 * Generate molecule from prompt
 */
export async function generate(prompt: string): Promise<GenerateResponse> {
  const response = await apiClient.post<GenerateResponse>('/generate', {
    prompt,
  })
  return response.data
}

/**
 * Fast prediction using ONNX model
 */
export async function predictFast(smiles: string): Promise<PredictResponse> {
  const response = await apiClient.post<PredictResponse>('/predict-fast', {
    smiles,
  })
  return response.data
}

/**
 * Molecule Library API
 */
export interface MoleculeItem {
  id: number
  name: string
  smiles?: string
  properties?: string
  thumbnail_b64?: string
  created_at: string
}

export interface MoleculeDetail extends MoleculeItem {
  json_graph?: string
  coords?: string
}

export interface SaveMoleculePayload {
  name: string
  smiles?: string
  json_graph?: string
  coords?: string
  properties?: string
  thumbnail_b64?: string | null
}

/**
 * Save molecule to library
 */
export async function saveMolecule(payload: SaveMoleculePayload): Promise<{ id: number; name: string; created_at: string }> {
  const response = await apiClient.post('/molecules/save', payload)
  return response.data
}

/**
 * List all molecules
 */
export async function listMolecules(limit: number = 50): Promise<MoleculeItem[]> {
  const response = await apiClient.get<MoleculeItem[]>(`/molecules/list?limit=${limit}`)
  return response.data
}

/**
 * Get molecule by ID
 */
export async function getMolecule(id: number): Promise<MoleculeDetail> {
  const response = await apiClient.get<MoleculeDetail>(`/molecules/${id}`)
  return response.data
}

/**
 * Delete molecule
 */
export async function deleteMolecule(id: number): Promise<{ status: string; id: number }> {
  const response = await apiClient.delete(`/molecules/${id}`)
  return response.data
}

