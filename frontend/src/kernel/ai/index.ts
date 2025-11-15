import axios from 'axios';

const API_URL = import.meta.env.VITE_API_URL || 'http://localhost:8000';

const apiClient = axios.create({
  baseURL: API_URL,
  headers: {
    'Content-Type': 'application/json',
  },
  timeout: 30000, // 30 seconds
});

interface RetryOptions {
  maxRetries?: number;
  retryDelay?: number;
  retryCondition?: (error: any) => boolean;
}

async function retryRequest<T>(
  request: () => Promise<T>,
  options: RetryOptions = {}
): Promise<T> {
  const {
    maxRetries = 3,
    retryDelay = 1000,
    retryCondition = (error) => {
      // Retry on network errors or 5xx errors
      if (!error.response) return true;
      return error.response.status >= 500;
    },
  } = options;
  
  let lastError: any;
  
  for (let attempt = 0; attempt <= maxRetries; attempt++) {
    try {
      return await request();
    } catch (error: any) {
      lastError = error;
      
      if (attempt < maxRetries && retryCondition(error)) {
        await new Promise((resolve) => setTimeout(resolve, retryDelay * (attempt + 1)));
        continue;
      }
      
      throw error;
    }
  }
  
  throw lastError;
}

export interface GenerateMoleculeProps {
  prompt: string;
  constraints?: {
    maxAtoms?: number;
    elements?: string[];
    properties?: {
      stability?: { min?: number; max?: number };
      toxicity?: { min?: number; max?: number };
    };
  };
}

export interface GenerateMoleculeResponse {
  smiles: string;
  properties?: {
    stability: number;
    toxicity: number;
    solubility: number;
    bioavailability: number;
    novelty: number;
  };
}

/**
 * Generate a molecule from a prompt
 */
export async function generateMolecule(
  props: GenerateMoleculeProps
): Promise<GenerateMoleculeResponse> {
  return retryRequest(async () => {
    const response = await apiClient.post<GenerateMoleculeResponse>('/api/v1/generate', {
      prompt: props.prompt,
      constraints: props.constraints,
    });
    return response.data;
  });
}

export interface MutateMoleculeProps {
  smiles: string;
  mutationType?: 'substitute' | 'add' | 'remove' | 'rearrange';
  target?: {
    atomIndex?: number;
    element?: string;
  };
}

export interface MutateMoleculeResponse {
  smiles: string;
  original_smiles: string;
  mutation_type: string;
  properties?: {
    stability: number;
    toxicity: number;
    solubility: number;
    bioavailability: number;
    novelty: number;
  };
}

/**
 * Mutate an existing molecule
 */
export async function mutateMolecule(
  props: MutateMoleculeProps
): Promise<MutateMoleculeResponse> {
  return retryRequest(async () => {
    const response = await apiClient.post<MutateMoleculeResponse>('/api/v1/mutate', {
      smiles: props.smiles,
      mutation_type: props.mutationType,
      target: props.target,
    });
    return response.data;
  });
}

export interface OptimizePropertiesProps {
  smiles: string;
  targetProperties: {
    stability?: number;
    toxicity?: number;
    solubility?: number;
    bioavailability?: number;
  };
  iterations?: number;
}

export interface OptimizePropertiesResponse {
  smiles: string;
  original_smiles: string;
  properties: {
    stability: number;
    toxicity: number;
    solubility: number;
    bioavailability: number;
    novelty: number;
  };
  optimization_history?: Array<{
    iteration: number;
    properties: Record<string, number>;
  }>;
}

/**
 * Optimize molecule properties
 */
export async function optimizeProperties(
  props: OptimizePropertiesProps
): Promise<OptimizePropertiesResponse> {
  return retryRequest(async () => {
    const response = await apiClient.post<OptimizePropertiesResponse>('/api/v1/optimize', {
      smiles: props.smiles,
      target_properties: props.targetProperties,
      iterations: props.iterations || 10,
    });
    return response.data;
  });
}

export interface ADMETProps {
  smiles: string;
}

export interface ADMETResponse {
  smiles: string;
  admet: {
    absorption: {
      score: number;
      caco2?: number;
      bioavailability?: number;
    };
    distribution: {
      score: number;
      vd?: number;
      ppb?: number;
    };
    metabolism: {
      score: number;
      cyp_substrates?: string[];
      half_life?: number;
    };
    excretion: {
      score: number;
      clearance?: number;
      renal_clearance?: number;
    };
    toxicity: {
      score: number;
      ld50?: number;
      mutagenic?: boolean;
      carcinogenic?: boolean;
    };
  };
}

/**
 * Predict ADMET properties
 */
export async function predictADMET(props: ADMETProps): Promise<ADMETResponse> {
  return retryRequest(async () => {
    const response = await apiClient.post<ADMETResponse>('/api/v1/admet', {
      smiles: props.smiles,
    });
    return response.data;
  });
}

