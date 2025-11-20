/**
 * User Molecules Store
 * Functions for accessing user_molecules table
 * Users can only access their own molecules (enforced by RLS)
 */
import { supabase, isSupabaseConfigured } from '../supabase';

export interface UserMolecule {
  id?: string;
  user_id: string;
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
 * Save a molecule to user_molecules for a specific user
 */
export async function saveUserMolecule(
  userId: string,
  molecule: Omit<UserMolecule, 'id' | 'user_id' | 'created_at' | 'updated_at'>
): Promise<string> {
  if (!isSupabaseConfigured() || !supabase) {
    throw new Error('Supabase is not configured. Please check your .env file and Supabase setup.');
  }
  
  const { data, error } = await supabase
    .from('user_molecules')
    .insert({
      ...molecule,
      user_id: userId,
      created_at: new Date().toISOString(),
      updated_at: new Date().toISOString(),
    })
    .select()
    .single();

  if (error) {
    throw new Error(`Failed to save user molecule: ${error.message}`);
  }

  return data.id;
}

/**
 * List all molecules for a specific user
 */
export async function listUserMolecules(userId: string): Promise<UserMolecule[]> {
  if (!isSupabaseConfigured() || !supabase) {
    throw new Error('Supabase is not configured. Please check your .env file and Supabase setup.');
  }
  
  const { data, error } = await supabase
    .from('user_molecules')
    .select('*')
    .eq('user_id', userId)
    .order('created_at', { ascending: false });

  if (error) {
    throw new Error(`Failed to list user molecules: ${error.message}`);
  }

  return data || [];
}

/**
 * Get a user molecule by ID
 */
export async function getUserMolecule(
  userId: string,
  moleculeId: string
): Promise<UserMolecule | null> {
  if (!isSupabaseConfigured() || !supabase) {
    throw new Error('Supabase is not configured. Please check your .env file and Supabase setup.');
  }
  
  const { data, error } = await supabase
    .from('user_molecules')
    .select('*')
    .eq('id', moleculeId)
    .eq('user_id', userId)
    .single();

  if (error) {
    if (error.code === 'PGRST116') {
      return null;
    }
    throw new Error(`Failed to get user molecule: ${error.message}`);
  }

  return data;
}

/**
 * Update a user molecule
 */
export async function updateUserMolecule(
  userId: string,
  moleculeId: string,
  updates: Partial<Omit<UserMolecule, 'id' | 'user_id' | 'created_at'>>
): Promise<void> {
  if (!isSupabaseConfigured() || !supabase) {
    throw new Error('Supabase is not configured. Please check your .env file and Supabase setup.');
  }
  
  const { error } = await supabase
    .from('user_molecules')
    .update({
      ...updates,
      updated_at: new Date().toISOString(),
    })
    .eq('id', moleculeId)
    .eq('user_id', userId);

  if (error) {
    throw new Error(`Failed to update user molecule: ${error.message}`);
  }
}

/**
 * Delete a user molecule
 */
export async function deleteUserMolecule(userId: string, moleculeId: string): Promise<void> {
  if (!isSupabaseConfigured() || !supabase) {
    throw new Error('Supabase is not configured. Please check your .env file and Supabase setup.');
  }
  
  const { error } = await supabase
    .from('user_molecules')
    .delete()
    .eq('id', moleculeId)
    .eq('user_id', userId);

  if (error) {
    throw new Error(`Failed to delete user molecule: ${error.message}`);
  }
}

/**
 * Search user molecules by name, formula, or SMILES
 */
export async function searchUserMolecules(
  userId: string,
  searchQuery: string
): Promise<UserMolecule[]> {
  if (!isSupabaseConfigured() || !supabase) {
    throw new Error('Supabase is not configured. Please check your .env file and Supabase setup.');
  }
  
  const query = searchQuery.toLowerCase().trim();
  if (!query) {
    return listUserMolecules(userId);
  }

  const { data, error } = await supabase
    .from('user_molecules')
    .select('*')
    .eq('user_id', userId)
    .or(`name.ilike.%${query}%,smiles.ilike.%${query}%,formula.ilike.%${query}%`)
    .order('created_at', { ascending: false });

  if (error) {
    throw new Error(`Failed to search user molecules: ${error.message}`);
  }

  return data || [];
}

/**
 * Fork a public molecule to user_molecules
 * Creates a copy that the user can edit
 */
export async function forkPublicMolecule(
  userId: string,
  publicMoleculeId: string
): Promise<string> {
  if (!isSupabaseConfigured() || !supabase) {
    throw new Error('Supabase is not configured. Please check your .env file and Supabase setup.');
  }
  
  // Get the public molecule
  const { data: publicMol, error: fetchError } = await supabase
    .from('public_molecules')
    .select('*')
    .eq('id', publicMoleculeId)
    .single();

  if (fetchError || !publicMol) {
    throw new Error('Failed to fetch public molecule for forking');
  }

  // Create a copy in user_molecules
  const { data, error } = await supabase
    .from('user_molecules')
    .insert({
      name: `${publicMol.name} (Forked)`,
      smiles: publicMol.smiles,
      formula: publicMol.formula,
      molfile: publicMol.molfile,
      thumbnail_b64: publicMol.thumbnail_b64,
      metadata: {
        ...publicMol.metadata,
        forked_from: publicMoleculeId,
        forked_at: new Date().toISOString(),
      },
      user_id: userId,
      created_at: new Date().toISOString(),
      updated_at: new Date().toISOString(),
    })
    .select()
    .single();

  if (error) {
    throw new Error(`Failed to fork molecule: ${error.message}`);
  }

  return data.id;
}

