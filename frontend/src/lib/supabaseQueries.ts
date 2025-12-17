/**
 * Supabase Queries with RLS Compliance
 * 
 * All queries assume:
 * - RLS is enabled on all tables
 * - User authentication is handled via Supabase Auth
 * - user_id comes from auth.uid() in RLS policies
 * 
 * These queries work with the anon key and respect RLS policies.
 */

import { supabase, isSupabaseConfigured } from '../supabase';

/**
 * FAVORITES TABLE QUERIES
 * RLS: Users can only SELECT/INSERT/DELETE their own favorites
 */

export interface Favorite {
  id: string;
  user_id: string;
  molecule_id: string;
  created_at: string;
}

/**
 * Add a favorite (requires authenticated user)
 * RLS policy ensures user_id matches auth.uid()
 */
export async function addFavorite(moleculeId: string): Promise<string> {
  if (!isSupabaseConfigured() || !supabase) {
    throw new Error('Supabase is not configured');
  }

  // Get current user - required for RLS
  const { data: { session } } = await supabase.auth.getSession();
  if (!session?.user) {
    throw new Error('Authentication required');
  }

  const { data, error } = await supabase
    .from('favorites')
    .insert({
      molecule_id: moleculeId,
      // user_id is automatically set by RLS policy (auth.uid())
      // But we can include it explicitly for clarity
    })
    .select()
    .single();

  if (error) {
    // Handle duplicate favorite case
    if (error.code === '23505') {
      throw new Error('Already favorited');
    }
    throw new Error(`Failed to add favorite: ${error.message}`);
  }

  return data.id;
}

/**
 * Remove a favorite (requires authenticated user)
 * RLS ensures user can only delete their own favorites
 */
export async function removeFavorite(moleculeId: string): Promise<void> {
  if (!isSupabaseConfigured() || !supabase) {
    throw new Error('Supabase is not configured');
  }

  const { data: { session } } = await supabase.auth.getSession();
  if (!session?.user) {
    throw new Error('Authentication required');
  }

  const { error } = await supabase
    .from('favorites')
    .delete()
    .eq('molecule_id', moleculeId);
    // RLS automatically filters by user_id = auth.uid()

  if (error) {
    throw new Error(`Failed to remove favorite: ${error.message}`);
  }
}

/**
 * List user's favorites (requires authenticated user)
 */
export async function listFavorites(): Promise<Favorite[]> {
  if (!isSupabaseConfigured() || !supabase) {
    throw new Error('Supabase is not configured');
  }

  const { data: { session } } = await supabase.auth.getSession();
  if (!session?.user) {
    throw new Error('Authentication required');
  }

  const { data, error } = await supabase
    .from('favorites')
    .select('*')
    .order('created_at', { ascending: false });
    // RLS automatically filters by user_id = auth.uid()

  if (error) {
    throw new Error(`Failed to list favorites: ${error.message}`);
  }

  return data || [];
}

/**
 * Check if a molecule is favorited by current user
 */
export async function isFavorited(moleculeId: string): Promise<boolean> {
  if (!isSupabaseConfigured() || !supabase) {
    return false;
  }

  const { data: { session } } = await supabase.auth.getSession();
  if (!session?.user) {
    return false;
  }

  const { data, error } = await supabase
    .from('favorites')
    .select('id')
    .eq('molecule_id', moleculeId)
    .maybeSingle();
    // RLS ensures we only check user's own favorites

  if (error && error.code !== 'PGRST116') {
    // PGRST116 = no rows returned (not an error for this use case)
    console.error('Error checking favorite:', error);
    return false;
  }

  return !!data;
}

/**
 * LAB_SESSIONS TABLE QUERIES
 * RLS: Users can only SELECT/INSERT/UPDATE/DELETE their own sessions
 */

export interface LabSession {
  id: string;
  user_id: string;
  name: string;
  scene_state: Record<string, any>;
  created_at: string;
  updated_at: string;
}

/**
 * Save a lab session (requires authenticated user)
 */
export async function saveLabSession(name: string, sceneState: Record<string, any>): Promise<string> {
  if (!isSupabaseConfigured() || !supabase) {
    throw new Error('Supabase is not configured');
  }

  const { data: { session } } = await supabase.auth.getSession();
  if (!session?.user) {
    throw new Error('Authentication required');
  }

  const { data, error } = await supabase
    .from('lab_sessions')
    .insert({
      name,
      scene_state: sceneState,
      // user_id handled by RLS
    })
    .select()
    .single();

  if (error) {
    throw new Error(`Failed to save lab session: ${error.message}`);
  }

  return data.id;
}

/**
 * Update a lab session (requires authenticated user)
 * RLS ensures user can only update their own sessions
 */
export async function updateLabSession(sessionId: string, sceneState: Record<string, any>): Promise<void> {
  if (!isSupabaseConfigured() || !supabase) {
    throw new Error('Supabase is not configured');
  }

  const { data: { session } } = await supabase.auth.getSession();
  if (!session?.user) {
    throw new Error('Authentication required');
  }

  const { error } = await supabase
    .from('lab_sessions')
    .update({
      scene_state: sceneState,
      updated_at: new Date().toISOString(),
    })
    .eq('id', sessionId);
    // RLS ensures user_id matches

  if (error) {
    throw new Error(`Failed to update lab session: ${error.message}`);
  }
}

/**
 * List user's lab sessions (requires authenticated user)
 */
export async function listLabSessions(): Promise<LabSession[]> {
  if (!isSupabaseConfigured() || !supabase) {
    throw new Error('Supabase is not configured');
  }

  const { data: { session } } = await supabase.auth.getSession();
  if (!session?.user) {
    throw new Error('Authentication required');
  }

  const { data, error } = await supabase
    .from('lab_sessions')
    .select('*')
    .order('updated_at', { ascending: false });
    // RLS automatically filters by user_id

  if (error) {
    throw new Error(`Failed to list lab sessions: ${error.message}`);
  }

  return data || [];
}

/**
 * USER_MOLECULES TABLE QUERIES
 * These already exist in userMoleculeStore.ts, but included here for reference
 * RLS: Users can only SELECT/INSERT/UPDATE/DELETE their own molecules
 * 
 * See userMoleculeStore.ts for implementation.
 */

/**
 * PUBLIC_MOLECULES TABLE
 * RLS: Public read access, no writes
 * 
 * These queries don't require authentication for reads.
 */

export interface PublicMolecule {
  id: string;
  name: string;
  smiles?: string;
  formula?: string;
  molfile?: string;
  thumbnail_b64?: string;
  created_at: string;
}

/**
 * List all public molecules (no auth required)
 */
export async function listPublicMolecules(limit: number = 50, offset: number = 0): Promise<PublicMolecule[]> {
  if (!isSupabaseConfigured() || !supabase) {
    throw new Error('Supabase is not configured');
  }

  const { data, error } = await supabase
    .from('public_molecules')
    .select('*')
    .order('created_at', { ascending: false })
    .range(offset, offset + limit - 1);

  if (error) {
    throw new Error(`Failed to list public molecules: ${error.message}`);
  }

  return data || [];
}

/**
 * Get a single public molecule by ID (no auth required)
 */
export async function getPublicMolecule(moleculeId: string): Promise<PublicMolecule | null> {
  if (!isSupabaseConfigured() || !supabase) {
    throw new Error('Supabase is not configured');
  }

  const { data, error } = await supabase
    .from('public_molecules')
    .select('*')
    .eq('id', moleculeId)
    .maybeSingle();

  if (error) {
    throw new Error(`Failed to get public molecule: ${error.message}`);
  }

  return data;
}

