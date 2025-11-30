import axios from 'axios'

const API_URL = import.meta.env.VITE_API_URL || 'http://localhost:8000'

export const apiClient = axios.create({
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
  molfile?: string
  formula?: string
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

/**
 * Convert SMILES to 3D molfile
 */
export async function convertSMILESToMolfile(smiles: string): Promise<{ molfile: string }> {
  const response = await apiClient.post<{ molfile: string }>('/convert/smiles', { smiles })
  return response.data
}

/**
 * Save molfile to molecule (persist after conversion)
 */
export async function saveMolfile(moleculeId: number, molfile: string): Promise<{ ok: boolean; id: number }> {
  const response = await apiClient.patch<{ ok: boolean; id: number }>(`/molecules/${moleculeId}/molfile`, { molfile })
  return response.data
}

/**
 * Thumbnail Generation API
 */
export interface GenerateThumbnailRequest {
  molfile?: string
  smiles?: string
  size?: [number, number]
  format?: 'PNG' | 'SVG'
}

export interface GenerateThumbnailBase64Request {
  molfile?: string
  smiles?: string
  size?: [number, number]
}

/**
 * Generate thumbnail image (returns image blob)
 */
export async function generateThumbnail(request: GenerateThumbnailRequest): Promise<Blob> {
  const response = await apiClient.post('/thumbnails/generate', request, {
    responseType: 'blob',
  })
  return response.data
}

/**
 * Generate thumbnail as base64 string
 */
export async function generateThumbnailBase64(
  request: GenerateThumbnailBase64Request
): Promise<{ thumbnail_b64: string }> {
  const response = await apiClient.post<{ thumbnail_b64: string }>(
    '/thumbnails/generate-base64',
    request
  )
  return response.data
}

/**
 * Generate optimized 3D molfile
 */
export async function generate3DMolfile(
  request: GenerateThumbnailBase64Request
): Promise<{ molfile: string }> {
  const response = await apiClient.post<{ molfile: string }>('/thumbnails/generate-3d', request)
  return response.data
}

/**
 * Admin Items API
 */
export interface Item {
  id: number
  name: string
  smiles?: string
  description?: string
  tags?: string[]
  status: 'in-stock' | 'sold-out' | 'archived'
  stock?: number
  structure_file_url?: string
  created_at: string
  updated_at: string
}

export interface CreateItemPayload {
  name: string
  smiles?: string
  description?: string
  tags?: string[]
  status?: 'in-stock' | 'sold-out' | 'archived'
  stock?: number
  structure_file?: File | null
}

export interface UpdateItemPayload extends Partial<CreateItemPayload> {}

/**
 * List all items
 */
export async function listItems(): Promise<Item[]> {
  const response = await apiClient.get<Item[]>('/api/v1/admin/items')
  return response.data
}

/**
 * Get item by ID
 */
export async function getItem(id: number): Promise<Item> {
  const response = await apiClient.get<Item>(`/api/v1/admin/items/${id}`)
  return response.data
}

/**
 * Create item
 */
export async function createItem(payload: CreateItemPayload): Promise<Item> {
  const formData = new FormData()
  formData.append('name', payload.name)
  if (payload.smiles) formData.append('smiles', payload.smiles)
  if (payload.description) formData.append('description', payload.description)
  if (payload.tags) formData.append('tags', JSON.stringify(payload.tags))
  if (payload.status) formData.append('status', payload.status)
  if (payload.stock !== undefined) formData.append('stock', payload.stock.toString())
  if (payload.structure_file) formData.append('structure_file', payload.structure_file)
  
  const response = await apiClient.post<Item>('/api/v1/admin/items', formData, {
    headers: {
      'Content-Type': 'multipart/form-data',
    },
  })
  return response.data
}

/**
 * Update item
 */
export async function updateItem(id: number, payload: UpdateItemPayload): Promise<Item> {
  const formData = new FormData()
  if (payload.name) formData.append('name', payload.name)
  if (payload.smiles !== undefined) formData.append('smiles', payload.smiles || '')
  if (payload.description !== undefined) formData.append('description', payload.description || '')
  if (payload.tags) formData.append('tags', JSON.stringify(payload.tags))
  if (payload.status) formData.append('status', payload.status)
  if (payload.stock !== undefined) formData.append('stock', payload.stock.toString())
  if (payload.structure_file) formData.append('structure_file', payload.structure_file)
  
  const response = await apiClient.put<Item>(`/api/v1/admin/items/${id}`, formData, {
    headers: {
      'Content-Type': 'multipart/form-data',
    },
  })
  return response.data
}

/**
 * Delete item
 */
export async function deleteItem(id: number): Promise<{ status: string; id: number }> {
  const response = await apiClient.delete(`/api/v1/admin/items/${id}`)
  return response.data
}

