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

