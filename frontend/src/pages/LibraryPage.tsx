/**
 * LibraryPage - Unified Library page with updated navigation
 * 
 * Features:
 * - Unified navigation tabs (All Molecules, My Library, Favorites)
 * - Favorites with star icons
 * - Sorting (A-Z, Z-A, Newest, Oldest)
 * - Search functionality
 * - Proper WebGL cleanup via unique keys
 */

import React, { useEffect, useMemo, useState, useRef } from 'react';
import { motion } from 'framer-motion';
import { supabase } from '../supabase';
import { useNavigate } from 'react-router-dom';
import {
  listUserMolecules,
  deleteUserMolecule,
  searchUserMolecules,
  forkPublicMolecule,
  updateUserMolecule,
  type UserMolecule
} from '../lib/userMoleculeStore';
import {
  listPublicMolecules,
  searchPublicMolecules,
  type PublicMolecule
} from '../lib/publicMoleculeStore';
import { convertSMILESToMolfile, saveMolfile } from '../lib/api';
import MoleculeCard from '../components/MoleculeCard';
import MoleculeImporter from '../components/MoleculeImporter';

type SortOption = 'a-z' | 'z-a' | 'newest' | 'oldest';
type ViewMode = 'all' | 'mine' | 'favorites';

export default function LibraryPage() {
  const navigate = useNavigate();
  const [activeView, setActiveView] = useState<ViewMode>('all');
  const [items, setItems] = useState<(PublicMolecule | UserMolecule)[]>([]);
  const [favorites, setFavorites] = useState<Set<string>>(new Set());
  const [loading, setLoading] = useState(false);
  const [q, setQ] = useState('');
  const [page, setPage] = useState(1);
  const [userId, setUserId] = useState<string | null>(null);
  const [sortBy, setSortBy] = useState<SortOption>('newest');
  const pageSize = 12;

  // Track which molecules have been converted to prevent duplicate conversions
  const convertedIdsRef = useRef<Set<string>>(new Set());
  const convertingRef = useRef<boolean>(false);

  // Load favorites from localStorage
  useEffect(() => {
    const stored = localStorage.getItem('molecule-favorites');
    if (stored) {
      try {
        const parsed = JSON.parse(stored);
        setFavorites(new Set(parsed));
      } catch (e) {
        console.warn('Failed to parse favorites:', e);
      }
    }
  }, []);

  // Save favorites to localStorage
  const toggleFavorite = (moleculeId: string) => {
    setFavorites((prev) => {
      const next = new Set(prev);
      if (next.has(moleculeId)) {
        next.delete(moleculeId);
      } else {
        next.add(moleculeId);
      }
      localStorage.setItem('molecule-favorites', JSON.stringify(Array.from(next)));
      return next;
    });
  };

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
    setLoading(true);
    setPage(1);

    try {
      if (activeView === 'mine') {
        // Load user molecules (auth required)
        if (!userId) {
          setItems([]);
          return;
        }
        const molecules = await listUserMolecules(userId);
        setItems(molecules);
      } else {
        // Load public molecules for 'all' and 'favorites' (default source)
        // Note: For favorites, we ideally want both sources, but for now we default to public
        // to maintain performance and simpler architecture.
        const molecules = await listPublicMolecules();
        setItems(molecules);
      }
    } catch (error) {
      console.error('Failed to load molecules:', error);
      setItems([]);
    } finally {
      setLoading(false);
    }
  }, [activeView, userId]);

  useEffect(() => {
    loadMolecules();

    // Set up realtime subscription for user_molecules table (only for 'mine' view)
    if (activeView === 'mine' && userId && supabase) {
      const channel = supabase
        .channel('user-molecules-changes')
        .on(
          'postgres_changes',
          {
            event: '*',
            schema: 'public',
            table: 'user_molecules',
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
  }, [activeView, userId, loadMolecules]);

  const openInLab = async (molecule: PublicMolecule | UserMolecule) => {
    try {
      const source = activeView === 'mine' ? 'user' : 'public';
      navigate(`/lab?id=${molecule.id}&source=${source}`, {
        state: {
          source,
          moleculeId: molecule.id,
          name: molecule.name,
          smiles: molecule.smiles,
          formula: (molecule as any).formula,
          molfile: (molecule as any).molfile,
        },
      });
    } catch (error) {
      console.error('Failed to open molecule:', error);
      alert('Failed to open molecule');
    }
  };

  const handleFork = async (molecule: PublicMolecule) => {
    if (!userId) {
      alert('Please sign in to fork molecules');
      return;
    }

    try {
      await forkPublicMolecule(userId, molecule.id!);
      alert(`"${molecule.name}" has been forked to your personal library!`);
      // If switching to 'mine' view context would be nice, but staying put is fine
      // User can switch to "My Library" to see it.
    } catch (error: any) {
      console.error('Failed to fork molecule:', error);
      alert(`Failed to fork molecule: ${error.message}`);
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
    const searchTimeout = setTimeout(async () => {
      if (q.trim()) {
        setLoading(true);
        try {
          if (activeView === 'mine') {
            if (!userId) {
              setItems([]);
              return;
            }
            const results = await searchUserMolecules(userId, q.trim());
            setItems(results);
          } else {
            // For 'all' and 'favorites'
            const results = await searchPublicMolecules(q.trim());
            setItems(results);
          }
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
  }, [q, activeView, userId, loadMolecules]);

  // Apply filtering and sorting
  const filtered = useMemo(() => {
    let result = [...items];

    // Apply favorites filter
    if (activeView === 'favorites') {
      result = result.filter((item) => item.id && favorites.has(String(item.id)));
    }

    // Apply sorting
    switch (sortBy) {
      case 'a-z':
        result.sort((a, b) => a.name.localeCompare(b.name));
        break;
      case 'z-a':
        result.sort((a, b) => b.name.localeCompare(a.name));
        break;
      case 'newest':
        result.sort((a, b) => new Date(b.created_at).getTime() - new Date(a.created_at).getTime());
        break;
      case 'oldest':
        result.sort((a, b) => new Date(a.created_at).getTime() - new Date(b.created_at).getTime());
        break;
    }

    return result;
  }, [items, activeView, favorites, sortBy]);

  const { data: paged, totalPages } = useMemo(() => {
    const start = (page - 1) * pageSize;
    const end = start + pageSize;
    const sliced = filtered.slice(start, end);
    const total = Math.ceil(filtered.length / pageSize);
    return {
      data: sliced,
      totalPages: total,
    };
  }, [filtered, page, pageSize]);

  // Reset page when filtered items change significantly
  useEffect(() => {
    const maxPage = Math.ceil(filtered.length / pageSize);
    if (maxPage > 0 && page > maxPage) {
      setPage(1);
    }
  }, [filtered.length, pageSize, page]);

  // Scroll to top when page changes
  useEffect(() => {
    window.scrollTo({ top: 0, behavior: 'smooth' });
  }, [page]);

  // Convert SMILES to molfile ONCE per molecule (runs after items are loaded)
  useEffect(() => {
    // Only convert if not already converting and items are loaded
    if (convertingRef.current || loading || items.length === 0) {
      return;
    }

    const convertMissingMolfiles = async () => {
      convertingRef.current = true;

      try {
        // Process molecules that need conversion
        const moleculesToConvert = items.filter(
          (item) =>
            !item.molfile &&
            item.smiles &&
            item.smiles.trim() &&
            item.id &&
            !convertedIdsRef.current.has(String(item.id))
        );

        if (moleculesToConvert.length === 0) {
          convertingRef.current = false;
          return;
        }

        console.log(`Converting ${moleculesToConvert.length} molecules...`);

        // Convert each molecule sequentially
        for (const molecule of moleculesToConvert) {
          const moleculeId = String(molecule.id);

          // Mark as being converted
          convertedIdsRef.current.add(moleculeId);

          try {
            // Convert SMILES to molfile
            const result = await convertSMILESToMolfile(molecule.smiles!);

            if (result.molfile && result.molfile.trim().length > 0) {
              // Update molecule in database (try both tables if generic, or specific based on view)
              // Since we know source from activeView, we can optimize:

              if (activeView === 'mine' && userId && typeof molecule.id === 'string') {
                await updateUserMolecule(userId, molecule.id, { molfile: result.molfile });
              } else if (typeof molecule.id === 'number') {
                // Numeric ID implies public table typically in this codebase?
                await saveMolfile(molecule.id, result.molfile);
              } else if (supabase && typeof molecule.id === 'string') {
                // Try public update (might fail RLS, but harmless)
                const { error } = await supabase
                  .from('public_molecules')
                  .update({ molfile: result.molfile })
                  .eq('id', molecule.id);
              }

              // Update local state WITHOUT triggering refetch
              setItems((prevItems) =>
                prevItems.map((item) =>
                  item.id === molecule.id
                    ? { ...item, molfile: result.molfile }
                    : item
                )
              );
            }
          } catch (error) {
            console.warn(`Failed to convert ${molecule.name}:`, error);
            convertedIdsRef.current.delete(moleculeId);
          }
        }
      } finally {
        convertingRef.current = false;
      }
    };

    const timeoutId = setTimeout(convertMissingMolfiles, 500);
    return () => clearTimeout(timeoutId);
  }, [items, activeView, userId, loading]);


  // Helper to determine if we are in 'mine' mode or 'public' mode for UI logic
  const isMine = activeView === 'mine';
  // Favorites view basically treats items as public-sourced for now
  const isPublicSource = !isMine;

  return (
    <motion.div
      initial={{ opacity: 0 }}
      animate={{ opacity: 1 }}
      transition={{ duration: 0.3 }}
      className="p-8 space-y-6"
    >
      <header className="mb-6">
        <div className="mb-4">
          <h1 className="text-3xl font-bold text-black truncate">Molecule Library</h1>
          <p className="text-darkGrey mt-1">
            Manage your molecules and explore the public collection
          </p>
        </div>

        {/* Removed old top buttons "Public Library" / "My Molecules" as requested */}

        {/* Search Bar with Refresh */}
        <div className="flex items-center gap-3 mb-4">
          <div className="flex-1 relative">
            <div className="absolute left-3 top-1/2 transform -translate-y-1/2 text-midGrey pointer-events-none">
              <svg width="20" height="20" viewBox="0 0 20 20" fill="none" xmlns="http://www.w3.org/2000/svg">
                <path d="M9 17C13.4183 17 17 13.4183 17 9C17 4.58172 13.4183 1 9 1C4.58172 1 1 4.58172 1 9C1 13.4183 4.58172 17 9 17Z" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" />
                <path d="M19 19L14.65 14.65" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" />
              </svg>
            </div>
            <input
              value={q}
              onChange={(e) => {
                setQ(e.target.value);
                setPage(1);
              }}
              placeholder={isMine ? "Search your library..." : "Search public molecules..."}
              className="w-full pl-10 pr-4 py-2 rounded-full border border-gray-200 shadow-sm focus:outline-none focus:ring-2 focus:ring-darkGrey/20 focus:border-darkGrey placeholder:text-midGrey"
            />
          </div>
          <button
            onClick={loadMolecules}
            className="p-2 rounded-full bg-white border border-gray-200 hover:bg-zinc-100 hover:border-zinc-300 active:bg-zinc-200 transition-all shadow-sm"
            title="Refresh"
          >
            <svg width="20" height="20" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
              <path d="M21.5 2V8M21.5 8H15.5M21.5 8L18 4.5C16.7429 3.24286 15.1767 2.34896 13.4606 1.91866C11.7446 1.48837 9.94453 1.53765 8.25383 2.06134C6.56313 2.58503 5.04987 3.56425 3.87691 4.89977C2.70395 6.23529 1.91462 7.87726 1.59 9.64M2.5 22V16M2.5 16H8.5M2.5 16L6 19.5C7.25714 20.7571 8.82332 21.651 10.5394 22.0813C12.2554 22.5116 14.0555 22.4624 15.7462 21.9387C17.4369 21.415 18.9501 20.4358 20.1231 19.1002C21.296 17.7647 22.0854 16.1227 22.41 14.36" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" />
            </svg>
          </button>
        </div>

        {/* Filter Tabs */}
        <div className="flex gap-2 mb-3 items-center justify-between">
          <div className="flex gap-2">
            <button
              className={`px-3 py-1.5 rounded-md text-sm font-medium transition-all ${activeView === 'all'
                ? 'bg-black text-white'
                : 'bg-zinc-100 text-zinc-700 hover:bg-zinc-200'
                }`}
              onClick={() => {
                setActiveView('all');
                setPage(1);
                setQ('');
              }}
            >
              All Molecules
            </button>
            <button
              className={`px-3 py-1.5 rounded-md text-sm font-medium transition-all ${activeView === 'mine'
                ? 'bg-black text-white'
                : 'bg-zinc-100 text-zinc-700 hover:bg-zinc-200'
                }`}
              onClick={() => {
                if (!userId) {
                  alert('Please sign in to view your library');
                  return;
                }
                setActiveView('mine');
                setPage(1);
                setQ('');
              }}
            >
              My Library
            </button>
            <button
              className={`px-3 py-1.5 rounded-md text-sm font-medium transition-all ${activeView === 'favorites'
                ? 'bg-black text-white'
                : 'bg-zinc-100 text-zinc-700 hover:bg-zinc-200'
                }`}
              onClick={() => {
                setActiveView('favorites');
                setPage(1);
                setQ('');
              }}
            >
              Favorites ({favorites.size})
            </button>
          </div>

          {/* Actions */}
          <div className="flex items-center gap-2">
            <MoleculeImporter userId={userId} onImportSuccess={loadMolecules} />
          </div>
        </div>

        {/* Sorting */}
        <div className="flex items-center gap-2">
          <span className="text-sm text-darkGrey">Sort by:</span>
          <select
            value={sortBy}
            onChange={(e) => setSortBy(e.target.value as SortOption)}
            className="px-3 py-1.5 rounded-md border border-gray-200 text-sm bg-white hover:border-gray-300 focus:outline-none focus:ring-2 focus:ring-blue-500/20"
          >
            <option value="newest">Newest → Oldest</option>
            <option value="oldest">Oldest → Newest</option>
            <option value="a-z">A → Z</option>
            <option value="z-a">Z → A</option>
          </select>
        </div>
      </header>

      {loading ? (
        <div className="text-center py-12 text-midGrey">Loading molecules...</div>
      ) : filtered.length === 0 ? (
        <div className="text-center py-12">
          <div className="text-midGrey mb-2">No results</div>
          <p className="text-midGrey text-sm">
            {q
              ? 'Try clearing the search'
              : activeView === 'favorites'
                ? 'No favorites yet. Click the star icon on molecule cards to add favorites.'
                : activeView === 'all'
                  ? 'No molecules in the public library yet'
                  : !userId
                    ? 'Please sign in to view your personal library'
                    : 'Save molecules in the Lab to see them here'}
          </p>
        </div>
      ) : (
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
          {paged.map((item) => {
            // For favorites view, we assume items are public for functionality purposes 
            // unless we add complexity.
            // If we are in 'mine' view, we show delete.
            // If we are in 'all' or 'favorites' view, we show fork.

            const moleculeId = String(item.id || `mol-${item.name}`);
            const isFavorite = favorites.has(moleculeId);

            return (
              <div key={`${activeView}-${moleculeId}`} className="relative">
                {/* Favorite Star */}
                <button
                  onClick={(e) => {
                    e.stopPropagation();
                    toggleFavorite(moleculeId);
                  }}
                  className="absolute top-2 right-2 z-20 p-2 rounded-full bg-white/90 hover:bg-white shadow-md transition-all hover:scale-110"
                  title={isFavorite ? 'Remove from favorites' : 'Add to favorites'}
                >
                  <svg
                    width="20"
                    height="20"
                    viewBox="0 0 24 24"
                    fill={isFavorite ? '#fbbf24' : 'none'}
                    stroke={isFavorite ? '#fbbf24' : '#6b7280'}
                    strokeWidth="2"
                    strokeLinecap="round"
                    strokeLinejoin="round"
                  >
                    <polygon points="12 2 15.09 8.26 22 9.27 17 14.14 18.18 21.02 12 17.77 5.82 21.02 7 14.14 2 9.27 8.91 8.26 12 2" />
                  </svg>
                </button>

                <MoleculeCard
                  item={{
                    id: moleculeId,
                    name: item.name,
                    smiles: item.smiles || undefined,
                    formula: item.formula || undefined,
                    molfile: item.molfile || undefined,
                    properties: (item as any).properties,
                    thumbnail_b64: item.thumbnail_b64 || undefined,
                    created_at: item.created_at,
                    user_id: isMine ? userId || undefined : undefined,
                  }}
                  onOpen={() => openInLab(item)}
                  onFork={isPublicSource ? () => handleFork(item as PublicMolecule) : undefined}
                  showFork={isPublicSource}
                  onDelete={isMine ? () => item.id && remove(String(item.id)) : undefined}
                />
              </div>
            );
          })}
        </div>
      )}

      {/* Pagination controls */}
      {filtered.length > 0 && totalPages > 1 && (
        <div className="flex items-center justify-center gap-2 pt-2">
          <button
            onClick={(e) => {
              e.preventDefault();
              e.stopPropagation();
              setPage((p) => Math.max(1, p - 1));
            }}
            disabled={page <= 1}
            className="btn-secondary px-3 py-1 text-sm disabled:opacity-50 disabled:cursor-not-allowed"
          >
            Prev
          </button>
          <div className="text-sm text-darkGrey">
            Page <span className="font-semibold text-black">{page}</span> of {totalPages}
          </div>
          <button
            onClick={(e) => {
              e.preventDefault();
              e.stopPropagation();
              if (page < totalPages) {
                setPage(page + 1);
              }
            }}
            disabled={page >= totalPages}
            className="btn-secondary px-3 py-1 text-sm disabled:opacity-50 disabled:cursor-not-allowed"
          >
            Next
          </button>
        </div>
      )}
    </motion.div>
  );
}
