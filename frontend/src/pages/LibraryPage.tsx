/**
 * LibraryPage - Unified Library page with Public and User tabs
 * 
 * Shows both public_molecules (read-only, forkable) and user_molecules (editable, deletable)
 */

import React, { useEffect, useMemo, useState, useCallback } from 'react';
import { motion } from 'framer-motion';
import { useInView } from 'react-intersection-observer';
import { supabase } from '../supabase';
import { 
  listUserMolecules, 
  deleteUserMolecule, 
  searchUserMolecules, 
  forkPublicMolecule,
  type UserMolecule 
} from '../lib/userMoleculeStore';
import { 
  listPublicMolecules, 
  searchPublicMolecules, 
  type PublicMolecule 
} from '../lib/publicMoleculeStore';
import MoleculeCard from '../components/MoleculeCard';

export default function LibraryPage() {
  const [tab, setTab] = useState<'public' | 'user'>('public');
  const [items, setItems] = useState<(PublicMolecule | UserMolecule)[]>([]);
  const [loading, setLoading] = useState(false);
  const [q, setQ] = useState('');
  const [page, setPage] = useState(1);
  const [userId, setUserId] = useState<string | null>(null);
  const [hasMore, setHasMore] = useState(true);
  const [isLoadingMore, setIsLoadingMore] = useState(false);
  const pageSize = 12;
  
  // Infinite scroll trigger
  const { ref: loadMoreRef, inView } = useInView({
    threshold: 0.1,
    rootMargin: '100px',
  });

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

  const loadMolecules = React.useCallback(async (reset: boolean = true) => {
    if (reset) {
      setLoading(true);
      setPage(1);
      setHasMore(true);
    } else {
      setIsLoadingMore(true);
    }
    
    try {
      if (tab === 'public') {
        // Load public molecules (no auth required)
        const molecules = await listPublicMolecules();
        if (reset) {
          setItems(molecules);
        } else {
          // For infinite scroll, append new items (though Supabase returns all)
          // In a real implementation with cursor-based pagination, you'd append here
          setItems(molecules);
        }
        setHasMore(false); // Supabase listPublicMolecules returns all, so no more pages
      } else {
        // Load user molecules (auth required)
        if (!userId) {
          setItems([]);
          setHasMore(false);
          return;
        }
        const molecules = await listUserMolecules(userId);
        if (reset) {
          setItems(molecules);
        } else {
          setItems(molecules);
        }
        setHasMore(false); // Supabase listUserMolecules returns all, so no more pages
      }
    } catch (error) {
      console.error('Failed to load molecules:', error);
      if (reset) {
        setItems([]);
      }
      setHasMore(false);
    } finally {
      if (reset) {
        setLoading(false);
      } else {
        setIsLoadingMore(false);
      }
    }
  }, [tab, userId]);

  useEffect(() => {
    loadMolecules();
    
    // Set up realtime subscription for user_molecules table (only for user tab)
    if (tab === 'user' && userId && supabase) {
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
  }, [tab, userId, loadMolecules]);

  const openInLab = async (molecule: PublicMolecule | UserMolecule) => {
    try {
      const source = tab === 'public' ? 'public' : 'user';
      window.location.href = `/lab?id=${molecule.id}&source=${source}`;
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
      // If on user tab, reload to show the forked molecule
      if (tab === 'user') {
        loadMolecules();
      }
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
          if (tab === 'public') {
            const results = await searchPublicMolecules(q.trim());
            setItems(results);
          } else {
            if (!userId) {
              setItems([]);
              return;
            }
            const results = await searchUserMolecules(userId, q.trim());
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
  }, [q, tab, userId, loadMolecules]);

  const filtered = useMemo(() => {
    // Items are already filtered by Supabase search or contain all molecules
    return items;
  }, [items]);

  const { data: paged, totalPages } = useMemo(() => {
    const start = (page - 1) * pageSize;
    const end = start + pageSize;
    const sliced = filtered.slice(start, end);
    const total = Math.ceil(filtered.length / pageSize);
    console.log('Pagination calc:', { 
      filteredLength: filtered.length, 
      page, 
      start, 
      end, 
      slicedLength: sliced.length,
      totalPages: total 
    });
    return {
      data: sliced,
      totalPages: total,
    };
  }, [filtered, page, pageSize]);

  // Reset page when filtered items change significantly (e.g., after search or tab change)
  useEffect(() => {
    const maxPage = Math.ceil(filtered.length / pageSize);
    if (maxPage > 0 && page > maxPage) {
      setPage(1);
    }
  }, [filtered.length, pageSize, page]);

  useEffect(() => {
    if (page > totalPages && totalPages > 0) {
      setPage(1);
    }
  }, [totalPages, page]);

  // Infinite scroll: load more when scroll trigger is in view
  useEffect(() => {
    if (inView && !loading && !isLoadingMore && page < totalPages) {
      const nextPage = page + 1;
      console.log('Infinite scroll triggered:', { current: page, next: nextPage, totalPages, filteredLength: filtered.length });
      setPage(nextPage);
    }
  }, [inView, loading, isLoadingMore, page, totalPages, filtered.length]);

  // Debug: Log molecule data
  useEffect(() => {
    if (paged.length > 0) {
      const sample = paged[0];
      const withMolfile = paged.filter(m => m.molfile && m.molfile.trim().length > 0).length;
      const withThumbnail = paged.filter(m => m.thumbnail_b64).length;
      
      console.log('Rendering molecules:', {
        total: paged.length,
        withMolfile,
        withThumbnail,
        sample: {
          name: sample.name,
          hasMolfile: !!sample.molfile,
          molfileLength: sample.molfile?.length || 0,
          molfilePreview: sample.molfile?.substring(0, 50) || 'null',
          hasThumbnail: !!sample.thumbnail_b64,
        }
      });
    }
  }, [paged]);

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
            {tab === 'public' 
              ? 'Browse curated molecules available to all users' 
              : 'Your personal saved molecular structures'}
          </p>
        </div>

        {/* Tabs */}
        <div className="flex gap-4 mb-4">
          <button
            className={`px-4 py-2 rounded-lg font-medium transition-all ${
              tab === 'public' 
                ? 'bg-blue-600 text-white shadow-md' 
                : 'bg-zinc-200 text-zinc-700 hover:bg-zinc-300'
            }`}
            onClick={() => {
              setTab('public');
              setPage(1);
              setQ('');
            }}
          >
            Public Library
          </button>
          <button
            className={`px-4 py-2 rounded-lg font-medium transition-all ${
              tab === 'user' 
                ? 'bg-blue-600 text-white shadow-md' 
                : 'bg-zinc-200 text-zinc-700 hover:bg-zinc-300'
            }`}
            onClick={() => {
              if (!userId) {
                alert('Please sign in to view your personal library');
                return;
              }
              setTab('user');
              setPage(1);
              setQ('');
            }}
          >
            My Molecules
          </button>
        </div>

        <div className="flex items-center gap-3">
          <div className="flex-1 relative">
            <div className="absolute left-3 top-1/2 transform -translate-y-1/2 text-midGrey pointer-events-none">
              <svg width="20" height="20" viewBox="0 0 20 20" fill="none" xmlns="http://www.w3.org/2000/svg">
                <path d="M9 17C13.4183 17 17 13.4183 17 9C17 4.58172 13.4183 1 9 1C4.58172 1 1 4.58172 1 9C1 13.4183 4.58172 17 9 17Z" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round"/>
                <path d="M19 19L14.65 14.65" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round"/>
              </svg>
            </div>
            <input
              value={q}
              onChange={(e) => {
                setQ(e.target.value);
                setPage(1);
              }}
              placeholder="Search by name, formula, or SMILES..."
              className="w-full pl-10 pr-4 py-2 rounded-full border border-gray-200 shadow-sm focus:outline-none focus:ring-2 focus:ring-darkGrey/20 focus:border-darkGrey placeholder:text-midGrey"
            />
          </div>
          <button
            onClick={loadMolecules}
            className="btn-secondary p-2 rounded-full hover:bg-zinc-200 transition-colors"
            title="Refresh"
          >
            <svg width="20" height="20" viewBox="0 0 20 20" fill="none" xmlns="http://www.w3.org/2000/svg">
              <path d="M17.5 2.5V7.5H12.5" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round"/>
              <path d="M2.5 10C2.5 13.5899 5.41015 16.5 9 16.5C10.8874 16.5 12.6082 15.7241 13.875 14.5L17.5 10.5" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round"/>
              <path d="M2.5 17.5L2.5 12.5H7.5" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round"/>
              <path d="M17.5 10C17.5 6.41015 14.5899 3.5 11 3.5C9.11258 3.5 7.39182 4.27588 6.125 5.5L2.5 9.5" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round"/>
            </svg>
          </button>
        </div>
      </header>

      {/* Info about 3D previews */}
      {filtered.length > 0 && (
        <div className="bg-blue-50 border border-blue-200 rounded-lg p-3 mb-4">
          <p className="text-sm text-blue-800">
            <strong>Note:</strong> 3D previews are always visible when molecules have molfile data. 
            Hover over a molecule to interact with it (rotate, zoom). If a molecule has SMILES but no molfile, 
            it will be automatically converted and saved. Thumbnails are shown as fallback.
          </p>
        </div>
      )}

      {loading ? (
        <div className="text-center py-12 text-midGrey">Loading molecules...</div>
      ) : filtered.length === 0 ? (
        <div className="text-center py-12">
          <div className="text-midGrey mb-2">No results</div>
          <p className="text-midGrey text-sm">
            {q 
              ? 'Try clearing the search' 
              : tab === 'public' 
                ? 'No molecules in the public library yet' 
                : !userId
                  ? 'Please sign in to view your personal library'
                  : 'Save molecules in the Lab to see them here'}
          </p>
        </div>
      ) : (
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
          {paged.map((item) => {
            const isPublic = tab === 'public';
            // Convert UUID string to number for MoleculeItem interface (use hash of UUID)
            const itemId = typeof item.id === 'string' 
              ? item.id.split('').reduce((acc, char) => acc + char.charCodeAt(0), 0) % 1000000
              : (item.id || 0);
            
            return (
              <MoleculeCard
                key={item.id || `mol-${item.name}`}
                item={{
                  id: typeof item.id === 'string' ? item.id : itemId,
                  name: item.name,
                  smiles: item.smiles || undefined,
                  formula: item.formula || undefined,
                  molfile: item.molfile || undefined,
                  properties: (item as any).properties,
                  thumbnail_b64: item.thumbnail_b64 || undefined,
                  created_at: item.created_at,
                  user_id: !isPublic && userId ? userId : undefined,
                }}
                onOpen={() => openInLab(item)}
                onFork={isPublic ? () => handleFork(item as PublicMolecule) : undefined}
                showFork={isPublic}
                onDelete={!isPublic && userId ? () => item.id && remove(item.id) : undefined}
                onMolfileUpdated={() => {
                  // Reload molecules to get updated molfile
                  loadMolecules();
                }}
              />
            );
          })}
        </div>
      )}

      {/* Infinite scroll trigger - must be after the grid */}
      {filtered.length > 0 && page < totalPages && (
        <div 
          ref={loadMoreRef} 
          className="h-20 w-full flex items-center justify-center py-4"
        >
          {isLoadingMore ? (
            <div className="text-sm text-midGrey flex items-center gap-2">
              <div className="w-4 h-4 border-2 border-midGrey border-t-transparent rounded-full animate-spin"></div>
              Loading more molecules...
            </div>
          ) : (
            <div className="text-xs text-midGrey">Scroll for more</div>
          )}
        </div>
      )}

      {/* Pagination controls (optional, can be hidden if using infinite scroll) */}
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
                const nextPage = page + 1;
                console.log('Next page clicked:', { current: page, next: nextPage, totalPages, filteredLength: filtered.length });
                setPage(nextPage);
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

