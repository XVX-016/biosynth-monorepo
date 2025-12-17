/**
 * Molecule Parser Utility
 * 
 * Converts molecules from Supabase format (molfile/json_graph) to Lab format.
 * Handles both molfile parsing and json_graph parsing.
 */

import { parseMolfile } from './molfileParser';
import { nanoid } from 'nanoid';
import type { Molecule, Atom, Bond } from '../types/molecule';
import { fromSmiles } from '../chemcore/io/smiles';

/**
 * Convert parsed molfile atoms/bonds to Lab format
 */
function convertMolfileToLab(
  parsed: { atoms: Array<{ x: number; y: number; z: number; element: string }>; bonds: Array<{ a1: number; a2: number; order?: number }> },
  metadata?: { id?: string; name?: string; smiles?: string; formula?: string; molfile?: string; source?: 'user' | 'public' }
): Molecule {
  const atomMap = new Map<number, string>(); // index -> atomId

  // Create atoms with IDs
  const atoms: Atom[] = parsed.atoms.map((atom, idx) => {
    const atomId = nanoid();
    atomMap.set(idx, atomId);
    // Validate position values
    if (typeof atom.x !== 'number' || typeof atom.y !== 'number' || typeof atom.z !== 'number') {
      console.warn('[convertMolfileToLab] Invalid atom position:', { atom, idx });
    }
    return {
      id: atomId,
      element: atom.element || 'C', // Default to carbon if missing
      position: { 
        x: typeof atom.x === 'number' ? atom.x : 0, 
        y: typeof atom.y === 'number' ? atom.y : 0, 
        z: typeof atom.z === 'number' ? atom.z : 0 
      },
    };
  });

  // Create bonds with IDs, referencing atom IDs
  const bonds: Bond[] = parsed.bonds
    .map((bond) => {
      const fromId = atomMap.get(bond.a1);
      const toId = atomMap.get(bond.a2);
      if (!fromId || !toId) return null;

      return {
        id: nanoid(),
        from: fromId,
        to: toId,
        order: (bond.order || 1) as 1 | 2 | 3,
      };
    })
    .filter((b): b is Bond => b !== null);

  return {
    id: metadata?.id || nanoid(),
    atoms,
    bonds,
    metadata: {
      name: metadata?.name,
      smiles: metadata?.smiles,
      formula: metadata?.formula,
      molfile: metadata?.molfile,
      source: metadata?.source,
    },
  };
}

/**
 * Convert json_graph to Lab format
 */
function convertJsonGraphToLab(
  jsonGraph: { atoms: Atom[]; bonds: Bond[] },
  metadata?: { id?: string; name?: string; smiles?: string; formula?: string; molfile?: string; source?: 'user' | 'public' }
): Molecule {
  // Ensure all atoms have IDs and proper position format
  const atoms: Atom[] = jsonGraph.atoms.map((atom, idx) => {
    const position = atom.position || { x: 0, y: 0, z: 0 };
    // Validate position format
    if (typeof position.x !== 'number' || typeof position.y !== 'number' || typeof position.z !== 'number') {
      console.warn('[convertJsonGraphToLab] Invalid atom position format:', { atom, idx });
      return {
        id: atom.id || nanoid(),
        element: atom.element || 'C',
        position: { x: 0, y: 0, z: 0 },
      };
    }
    return {
      id: atom.id || nanoid(),
      element: atom.element || 'C',
      position,
    };
  });

  // Ensure all bonds have IDs and proper order
  const bonds: Bond[] = jsonGraph.bonds.map((bond, idx) => {
    // Validate bond references
    const fromAtom = atoms.find(a => a.id === bond.from);
    const toAtom = atoms.find(a => a.id === bond.to);
    if (!fromAtom || !toAtom) {
      console.warn('[convertJsonGraphToLab] Bond references invalid atom:', { bond, idx, fromAtom: !!fromAtom, toAtom: !!toAtom });
      return null;
    }
    return {
      id: bond.id || nanoid(),
      from: bond.from,
      to: bond.to,
      order: (bond.order || 1) as 1 | 2 | 3,
    };
  }).filter((b): b is Bond => b !== null);

  // Validate molecule structure
  if (atoms.length === 0 && bonds.length > 0) {
    console.warn('[convertJsonGraphToLab] Invalid: bonds without atoms');
  }

  return {
    id: metadata?.id || nanoid(),
    atoms,
    bonds,
    metadata: {
      name: metadata?.name,
      smiles: metadata?.smiles,
      formula: metadata?.formula,
      molfile: metadata?.molfile,
      source: metadata?.source,
    },
  };
}

/**
 * Convert SMILES to Lab format using frontend parser
 */
function convertSmilesToLab(
  smiles: string,
  metadata?: { id?: string; name?: string; smiles?: string; formula?: string; molfile?: string; source?: 'user' | 'public' }
): Molecule {
  try {
    const parsed = fromSmiles(smiles);
    // Ensure positions have valid z values (OpenChemLib may return 2D with z=0)
    const atoms: Atom[] = parsed.atoms.map((atom) => ({
      ...atom,
      position: {
        x: atom.position?.x || 0,
        y: atom.position?.y || 0,
        z: atom.position?.z || 0, // Default to 0 if missing (2D layout)
      },
    }));

    return {
      id: metadata?.id || nanoid(),
      atoms,
      bonds: parsed.bonds,
      metadata: {
        name: metadata?.name,
        smiles: metadata?.smiles || smiles,
        formula: metadata?.formula,
        molfile: metadata?.molfile,
        source: metadata?.source,
      },
    };
  } catch (e) {
    console.warn('[convertSmilesToLab] Failed to parse SMILES:', e);
    throw e;
  }
}

/**
 * Parse molecule from Supabase data
 * 
 * Priority:
 * 1. json_graph (if available and valid)
 * 2. molfile (if available)
 * 3. SMILES (if available, generate structure)
 * 4. Return empty molecule with metadata
 */
export function parseMoleculeFromSupabase(data: {
  id?: string | number;
  name?: string;
  smiles?: string | null;
  formula?: string | null;
  molfile?: string | null;
  json_graph?: string | null;
  source?: 'user' | 'public';
}): Molecule {
  const metadata = {
    id: String(data.id || nanoid()),
    name: data.name,
    smiles: data.smiles || undefined,
    formula: data.formula || undefined,
    molfile: data.molfile || undefined,
    source: data.source,
  };

  // Try json_graph first
  if (data.json_graph) {
    try {
      const parsed = typeof data.json_graph === 'string' 
        ? JSON.parse(data.json_graph) 
        : data.json_graph;
      
      if (parsed && typeof parsed === 'object' && 'atoms' in parsed && 'bonds' in parsed) {
        const molecule = convertJsonGraphToLab(parsed, metadata);
        // Validate parsed result
        if (molecule.atoms.length > 0) {
          console.log('[parseMoleculeFromSupabase] Successfully parsed json_graph:', {
            atomCount: molecule.atoms.length,
            bondCount: molecule.bonds.length,
          });
          return molecule;
        } else {
          console.warn('[parseMoleculeFromSupabase] json_graph parsed but resulted in empty atoms');
        }
      } else {
        console.warn('[parseMoleculeFromSupabase] json_graph missing required fields (atoms/bonds)');
      }
    } catch (e) {
      console.warn('[parseMoleculeFromSupabase] Failed to parse json_graph:', e);
    }
  }

  // Fallback to molfile
  if (data.molfile) {
    try {
      const parsed = parseMolfile(data.molfile);
      if (parsed.atoms.length > 0) {
        const molecule = convertMolfileToLab(parsed, metadata);
        console.log('[parseMoleculeFromSupabase] Successfully parsed molfile:', {
          atomCount: molecule.atoms.length,
          bondCount: molecule.bonds.length,
        });
        return molecule;
      } else {
        console.warn('[parseMoleculeFromSupabase] molfile parsed but resulted in empty atoms');
      }
    } catch (e) {
      console.warn('[parseMoleculeFromSupabase] Failed to parse molfile:', e);
    }
  }

  // Fallback to SMILES (generate structure from SMILES)
  if (data.smiles && data.smiles.trim()) {
    try {
      const molecule = convertSmilesToLab(data.smiles.trim(), metadata);
      if (molecule.atoms.length > 0) {
        console.log('[parseMoleculeFromSupabase] Generated from SMILES:', {
          atomCount: molecule.atoms.length,
          bondCount: molecule.bonds.length,
          name: metadata.name,
        });
        return molecule;
      } else {
        console.warn('[parseMoleculeFromSupabase] SMILES parsed but resulted in empty atoms');
      }
    } catch (e) {
      console.warn('[parseMoleculeFromSupabase] Failed to parse SMILES:', e);
    }
  }

  // Return empty molecule with metadata
  const errorMsg = data.smiles
    ? `No structural data (molfile/json_graph) and SMILES conversion failed for "${metadata.name || 'unknown molecule'}"`
    : `No valid molecule data found for "${metadata.name || 'unknown molecule'}" (missing molfile, json_graph, and SMILES)`;
  console.warn(`[parseMoleculeFromSupabase] ${errorMsg}`);
  return {
    id: metadata.id,
    atoms: [],
    bonds: [],
    metadata,
  };
}

