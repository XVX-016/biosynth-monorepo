import { supabase, isSupabaseConfigured } from '../supabase';

export interface SupabaseMolecule {
  id?: string;
  name: string;
  smiles?: string;
  formula?: string;
  molfile?: string; // V2000/V3000 molfile for 3D coordinates
  json_graph?: string;
  properties?: string;
  thumbnail_b64?: string;
  user_id: string;
  created_at: string;
}

/**
 * Save a molecule to Supabase for a specific user
 */
export async function saveMolecule(userId: string, molecule: Omit<SupabaseMolecule, 'id' | 'user_id' | 'created_at'>): Promise<string> {
  if (!isSupabaseConfigured() || !supabase) {
    throw new Error('Supabase is not configured. Please check your .env file and Supabase setup.');
  }

  const { data, error } = await supabase
    .from('molecules')
    .insert({
      ...molecule,
      user_id: userId,
      created_at: new Date().toISOString(),
    })
    .select()
    .single();

  if (error) {
    throw new Error(`Failed to save molecule: ${error.message}`);
  }

  return data.id;
}

/**
 * List all molecules for a specific user
 */
export async function listMolecules(userId: string): Promise<SupabaseMolecule[]> {
  if (!isSupabaseConfigured() || !supabase) {
    throw new Error('Supabase is not configured. Please check your .env file and Supabase setup.');
  }

  const { data, error } = await supabase
    .from('molecules')
    .select('*')
    .eq('user_id', userId)
    .order('created_at', { ascending: false });

  if (error) {
    throw new Error(`Failed to list molecules: ${error.message}`);
  }

  return data || [];
}

/**
 * List ALL molecules without user filter
 */
export async function listAllMolecules(): Promise<SupabaseMolecule[]> {
  if (!isSupabaseConfigured() || !supabase) {
    console.warn('Supabase is not configured. Please check your .env file and Supabase setup.');
    return [];
  }

  const { data, error } = await supabase
    .from('molecules')
    .select('*')
    .order('created_at', { ascending: false });

  if (error) {
    throw new Error(`Failed to list all molecules: ${error.message}`);
  }

  return data || [];
}

/**
 * Delete a molecule for a specific user
 */
export async function deleteMolecule(userId: string, moleculeId: string): Promise<void> {
  if (!isSupabaseConfigured() || !supabase) {
    throw new Error('Supabase is not configured. Please check your .env file and Supabase setup.');
  }

  const { error } = await supabase
    .from('molecules')
    .delete()
    .eq('id', moleculeId)
    .eq('user_id', userId);

  if (error) {
    throw new Error(`Failed to delete molecule: ${error.message}`);
  }
}

/**
 * Get a molecule by ID for a specific user
 */
export async function getMolecule(userId: string, moleculeId: string): Promise<SupabaseMolecule | null> {
  if (!isSupabaseConfigured() || !supabase) {
    throw new Error('Supabase is not configured. Please check your .env file and Supabase setup.');
  }

  const { data, error } = await supabase
    .from('molecules')
    .select('*')
    .eq('id', moleculeId)
    .eq('user_id', userId)
    .single();

  if (error) {
    if (error.code === 'PGRST116') {
      // No rows returned
      return null;
    }
    throw new Error(`Failed to get molecule: ${error.message}`);
  }

  return data;
}

/**
 * Search molecules by name or SMILES for a specific user
 */
export async function searchMolecules(
  userId: string,
  searchQuery: string
): Promise<SupabaseMolecule[]> {
  if (!isSupabaseConfigured() || !supabase) {
    throw new Error('Supabase is not configured. Please check your .env file and Supabase setup.');
  }

  const query = searchQuery.toLowerCase().trim();
  if (!query) {
    return listMolecules(userId);
  }

  const { data, error } = await supabase
    .from('molecules')
    .select('*')
    .eq('user_id', userId)
    .or(`name.ilike.%${query}%,smiles.ilike.%${query}%,formula.ilike.%${query}%`)
    .order('created_at', { ascending: false });

  if (error) {
    throw new Error(`Failed to search molecules: ${error.message}`);
  }

  return data || [];
}

/**
 * Update a molecule for a specific user
 */
export async function updateMolecule(
  userId: string,
  moleculeId: string,
  updates: Partial<Omit<SupabaseMolecule, 'id' | 'user_id' | 'created_at'>>
): Promise<void> {
  if (!isSupabaseConfigured() || !supabase) {
    throw new Error('Supabase is not configured. Please check your .env file and Supabase setup.');
  }

  const { error } = await supabase
    .from('molecules')
    .update(updates)
    .eq('id', moleculeId)
    .eq('user_id', userId);

  if (error) {
    throw new Error(`Failed to update molecule: ${error.message}`);
  }
}

