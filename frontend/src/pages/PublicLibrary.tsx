/**
 * PublicLibrary - Displays the global public_molecules library
 * 
 * Shows read-only molecules available to all users
 * Supports search and filtering
 */
import React, { useEffect, useState, useMemo } from 'react';
import { motion } from 'framer-motion';
import { listPublicMolecules, searchPublicMolecules, filterPublicMoleculesByElement, type PublicMolecule } from '../lib/publicMoleculeStore';
import MoleculeCard from '../components/MoleculeCard';
import MoleculeFilters, { type MoleculeFilters as Filters } from '../components/MoleculeFilters';
import { forkPublicMolecule } from '../lib/userMoleculeStore';
import { supabase } from '../supabase';
import { useNavigate } from 'react-router-dom';

export default function PublicLibrary() {
  const navigate = useNavigate();
  const [molecules, setMolecules] = useState<PublicMolecule[]>([]);
  const [loading, setLoading] = useState(false);
  const [filters, setFilters] = useState<Filters>({ query: '', element: '' });
  const [page, setPage] = useState(1);
  const [userId, setUserId] = useState<string | null>(null);
  const pageSize = 12;

  useEffect(() => {
    if (!supabase) return;

    supabase.auth.getSession().then(({ data: { session } }) => {
      if (session?.user) {
        setUserId(session.user.id);
      }
    });

    const {
      data: { subscription },
    } = supabase.auth.onAuthStateChange((_event, session) => {
      if (session?.user) {
        setUserId(session.user.id);
      } else {
        setUserId(null);
      }
    });

    return () => subscription.unsubscribe();
  }, []);

  const loadMolecules = React.useCallback(async () => {
    setLoading(true);
    try {
      let results: PublicMolecule[];
      
      if (filters.query) {
        results = await searchPublicMolecules(filters.query);
      } else if (filters.element) {
        results = await filterPublicMoleculesByElement(filters.element);
      } else {
        results = await listPublicMolecules();
      }

      // Apply element filter if both query and element are set
      if (filters.query && filters.element) {
        results = results.filter(mol => 
          mol.formula?.toUpperCase().includes(filters.element.toUpperCase())
        );
      }

      setMolecules(results);
    } catch (error) {
      console.error('Failed to load public molecules:', error);
    } finally {
      setLoading(false);
    }
  }, [filters]);

  useEffect(() => {
    loadMolecules();
  }, [loadMolecules]);

  const handleFork = async (molecule: PublicMolecule) => {
    if (!userId) {
      alert('Please sign in to fork molecules');
      return;
    }

    try {
      await forkPublicMolecule(userId, molecule.id!);
      // Show toast notification (you can replace with a proper toast library)
      alert(`"${molecule.name}" has been forked to your personal library!`);
    } catch (error: any) {
      console.error('Failed to fork molecule:', error);
      alert(`Failed to fork molecule: ${error.message}`);
    }
  };

  const openInLab = async (molecule: PublicMolecule) => {
    try {
      const source = 'public';
      navigate(`/lab?id=${molecule.id}&source=${source}`, {
        state: {
          source,
          moleculeId: molecule.id,
          name: molecule.name,
          smiles: molecule.smiles,
          formula: molecule.formula,
          molfile: molecule.molfile,
        },
      });
    } catch (error) {
      console.error('Failed to open molecule:', error);
      alert('Failed to open molecule');
    }
  };

  const { paged, totalPages } = useMemo(() => {
    const start = (page - 1) * pageSize;
    const end = start + pageSize;
    return {
      paged: molecules.slice(start, end),
      totalPages: Math.ceil(molecules.length / pageSize),
    };
  }, [molecules, page, pageSize]);

  useEffect(() => {
    if (page > totalPages) setPage(1);
  }, [totalPages, page]);

  return (
    <motion.div
      initial={{ opacity: 0 }}
      animate={{ opacity: 1 }}
      transition={{ duration: 0.3 }}
      className="p-8 space-y-6"
    >
      <header>
        <h1 className="text-3xl font-bold text-black mb-2">Global Molecule Library</h1>
        <p className="text-darkGrey">Browse curated molecules available to all users</p>
      </header>

      <div className="mb-6">
        <MoleculeFilters value={filters} onChange={setFilters} />
      </div>

      {loading ? (
        <div className="text-center py-12 text-midGrey">Loading molecules...</div>
      ) : molecules.length === 0 ? (
        <div className="text-center py-12">
          <div className="text-midGrey mb-2">No molecules found</div>
          <p className="text-midGrey text-sm">
            {filters.query || filters.element ? 'Try adjusting your search filters' : 'No molecules in the public library yet'}
          </p>
        </div>
      ) : (
        <>
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
            {paged.map((molecule) => (
              <MoleculeCard
                key={molecule.id}
                item={{
                  id: parseInt(molecule.id || '0'),
                  name: molecule.name,
                  smiles: molecule.smiles,
                  formula: molecule.formula,
                  molfile: molecule.molfile,
                  thumbnail_b64: molecule.thumbnail_b64,
                  created_at: molecule.created_at,
                }}
                onOpen={() => openInLab(molecule)}
                onFork={() => handleFork(molecule)}
                showFork={true}
                onDelete={undefined} // Public molecules can't be deleted
              />
            ))}
          </div>

          {/* Pagination */}
          {molecules.length > 0 && totalPages > 1 && (
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
        </>
      )}
    </motion.div>
  );
}

