/**
 * Public Molecules Store
 * Functions for accessing the global public_molecules table
 * These molecules are read-only for regular users, editable only by admins
 */
import { supabase, isSupabaseConfigured } from '../supabase';

export interface PublicMolecule {
  id?: string;
  name: string;
  smiles?: string;
  formula?: string;
  molfile?: string;
  thumbnail_b64?: string;
  metadata?: Record<string, any>;
  created_at: string;
  updated_at?: string;
}

/**
 * List all public molecules (read-only for all users)
 */
export async function listPublicMolecules(): Promise<PublicMolecule[]> {
  if (!isSupabaseConfigured() || !supabase) {
    throw new Error('Supabase is not configured. Please check your .env file and Supabase setup.');
  }
  
  const { data, error } = await supabase
    .from('public_molecules')
    .select('*')
    .order('name', { ascending: true });

  if (error) {
    throw new Error(`Failed to list public molecules: ${error.message}`);
  }

  return data || [];
}

/**
 * Get a public molecule by ID
 */
export async function getPublicMolecule(moleculeId: string): Promise<PublicMolecule | null> {
  if (!isSupabaseConfigured() || !supabase) {
    throw new Error('Supabase is not configured. Please check your .env file and Supabase setup.');
  }
  
  const { data, error } = await supabase
    .from('public_molecules')
    .select('*')
    .eq('id', moleculeId)
    .single();

  if (error) {
    if (error.code === 'PGRST116') {
      return null;
    }
    throw new Error(`Failed to get public molecule: ${error.message}`);
  }

  return data;
}

/**
 * Search public molecules by name, formula, or SMILES
 */
export async function searchPublicMolecules(searchQuery: string): Promise<PublicMolecule[]> {
  if (!isSupabaseConfigured() || !supabase) {
    throw new Error('Supabase is not configured. Please check your .env file and Supabase setup.');
  }
  
  const query = searchQuery.toLowerCase().trim();
  if (!query) {
    return listPublicMolecules();
  }

  const { data, error } = await supabase
    .from('public_molecules')
    .select('*')
    .or(`name.ilike.%${query}%,smiles.ilike.%${query}%,formula.ilike.%${query}%`)
    .order('name', { ascending: true });

  if (error) {
    throw new Error(`Failed to search public molecules: ${error.message}`);
  }

  return data || [];
}

/**
 * Filter public molecules by element (checks formula)
 */
export async function filterPublicMoleculesByElement(element: string): Promise<PublicMolecule[]> {
  if (!isSupabaseConfigured() || !supabase) {
    throw new Error('Supabase is not configured. Please check your .env file and Supabase setup.');
  }
  
  if (!element || element.trim() === '') {
    return listPublicMolecules();
  }

  const { data, error } = await supabase
    .from('public_molecules')
    .select('*')
    .ilike('formula', `%${element}%`)
    .order('name', { ascending: true });

  if (error) {
    throw new Error(`Failed to filter public molecules: ${error.message}`);
  }

  return data || [];
}

