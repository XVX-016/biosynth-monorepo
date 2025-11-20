/**
 * LibraryPage - Alternative Library page implementation using Supabase
 * 
 * This is an optional implementation that uses Supabase instead of the backend API.
 * To use this, update your routing to point to LibraryPage instead of Library.
 */

import React, { useEffect, useMemo, useState } from 'react';
import { motion } from 'framer-motion';
import { supabase } from '../supabase';
import { listUserMolecules, deleteUserMolecule, searchUserMolecules, type UserMolecule } from '../lib/userMoleculeStore';
import MoleculeCard from '../components/MoleculeCard';
import SearchBar from '../components/SearchBar';

export default function LibraryPage() {
  const [items, setItems] = useState<UserMolecule[]>([]);
  const [loading, setLoading] = useState(false);
  const [q, setQ] = useState('');
  const [page, setPage] = useState(1);
  const [userId, setUserId] = useState<string | null>(null);
  const pageSize = 12;

  useEffect(() => {
    if (!supabase) return;

    // Check current session
    supabase.auth.getSession().then(({ data: { session } }) => {
      if (session?.user) {
        setUserId(session.user.id);
      } else {
        setUserId(null);
        setItems([]);
      }
    });

    // Listen for auth changes
    const {
      data: { subscription },
    } = supabase.auth.onAuthStateChange((_event, session) => {
      if (session?.user) {
        setUserId(session.user.id);
      } else {
        setUserId(null);
        setItems([]);
      }
    });

    return () => subscription.unsubscribe();
  }, []);

  const loadMolecules = React.useCallback(async () => {
    if (!userId) return;
    
    setLoading(true);
    try {
      const molecules = await listUserMolecules(userId);
      setItems(molecules);
    } catch (error) {
      console.error('Failed to load molecules:', error);
    } finally {
      setLoading(false);
    }
  }, [userId]);

  useEffect(() => {
    if (userId) {
      loadMolecules();
      
      // Set up realtime subscription for molecules table
      if (supabase) {
        const channel = supabase
          .channel('molecules-changes')
          .on(
            'postgres_changes',
            {
              event: '*',
              schema: 'public',
              table: 'molecules',
              filter: `user_id=eq.${userId}`,
            },
            () => {
              // Reload molecules when changes occur
              loadMolecules();
            }
          )
          .subscribe();

        return () => {
          supabase.removeChannel(channel);
        };
      }
    }
  }, [userId, loadMolecules]);

  const openInLab = async (molecule: UserMolecule) => {
    try {
      // Navigate to lab with molecule ID and source
      window.location.href = `/lab?id=${molecule.id}&source=user`;
    } catch (error) {
      console.error('Failed to open molecule:', error);
      alert('Failed to open molecule');
    }
  };

  const remove = async (moleculeId: string) => {
    if (!userId) return;
    if (!confirm('Delete this molecule?')) return;
    
    try {
      await deleteUserMolecule(userId, moleculeId);
      setItems(items.filter((i) => i.id !== moleculeId));
    } catch (error) {
      console.error('Failed to delete molecule:', error);
      alert('Failed to delete molecule');
    }
  };

  // Use Supabase search when query is provided
  useEffect(() => {
    if (!userId) return;
    
    const searchTimeout = setTimeout(async () => {
      if (q.trim()) {
        setLoading(true);
        try {
          const results = await searchUserMolecules(userId, q.trim());
          setItems(results);
        } catch (error) {
          console.error('Search error:', error);
          // Fallback: reload all molecules on error
          loadMolecules();
        } finally {
          setLoading(false);
        }
      } else {
        // If search is cleared, reload all molecules
        loadMolecules();
      }
    }, 300); // Debounce search

    return () => clearTimeout(searchTimeout);
  }, [q, userId, loadMolecules]);

  const filtered = useMemo(() => {
    // Items are already filtered by Supabase search or contain all molecules
    return items;
  }, [items]);

  const { data: paged, totalPages } = useMemo(() => {
    const start = (page - 1) * pageSize;
    const end = start + pageSize;
    return {
      data: filtered.slice(start, end),
      totalPages: Math.ceil(filtered.length / pageSize),
    };
  }, [filtered, page]);

  useEffect(() => {
    if (page > totalPages) setPage(1);
  }, [totalPages, page]);

  if (!userId) {
    return (
      <motion.div
        initial={{ opacity: 0 }}
        animate={{ opacity: 1 }}
        transition={{ duration: 0.3 }}
        className="p-8 space-y-6"
      >
        <div className="text-center py-12">
          <h2 className="text-2xl font-bold text-black mb-2">Authentication Required</h2>
          <p className="text-midGrey">Please sign in to view your molecule library.</p>
        </div>
      </motion.div>
    );
  }

  return (
    <motion.div
      initial={{ opacity: 0 }}
      animate={{ opacity: 1 }}
      transition={{ duration: 0.3 }}
      className="p-8 space-y-6"
    >
      <header className="mb-6">
        <div className="mb-4">
          <h1 className="text-3xl font-bold text-black truncate">My Molecule Library</h1>
          <p className="text-darkGrey mt-1">Your personal saved molecular structures</p>
        </div>
        <div className="flex items-center gap-3">
          <SearchBar
            value={q}
            onChange={(query) => {
              setQ(query);
              setPage(1);
            }}
            placeholder="Search by name, formula, or SMILES..."
          />
          <button
            onClick={loadMolecules}
            className="btn-secondary px-4 py-2 whitespace-nowrap"
          >
            Refresh
          </button>
        </div>
      </header>

      {loading ? (
        <div className="text-center py-12 text-midGrey">Loading molecules...</div>
      ) : filtered.length === 0 ? (
        <div className="text-center py-12">
          <div className="text-midGrey mb-2">No results</div>
          <p className="text-midGrey text-sm">
            {q ? 'Try clearing the search' : 'Save molecules in the Lab to see them here'}
          </p>
        </div>
      ) : (
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 xl:grid-cols-4 gap-6">
          {paged.map((item) => (
            <MoleculeCard
              key={item.id}
              item={{
                id: parseInt(item.id || '0'),
                name: item.name,
                smiles: item.smiles,
                formula: item.formula,
                molfile: item.molfile,
                properties: item.properties,
                thumbnail_b64: item.thumbnail_b64,
                created_at: item.created_at,
              }}
              onOpen={() => openInLab(item)}
              onDelete={() => item.id && remove(item.id)}
              showFork={false}
            />
          ))}
        </div>
      )}

      {/* Pagination */}
      {filtered.length > 0 && totalPages > 1 && (
        <div className="flex items-center justify-center gap-2 pt-2">
          <button
            onClick={() => setPage((p) => Math.max(1, p - 1))}
            disabled={page === 1}
            className="btn-secondary px-3 py-1 text-sm disabled:opacity-50"
          >
            Prev
          </button>
          <div className="text-sm text-darkGrey">
            Page <span className="font-semibold text-black">{page}</span> of {totalPages}
          </div>
          <button
            onClick={() => setPage((p) => Math.min(totalPages, p + 1))}
            disabled={page === totalPages}
            className="btn-secondary px-3 py-1 text-sm disabled:opacity-50"
          >
            Next
          </button>
        </div>
      )}
    </motion.div>
  );
}

