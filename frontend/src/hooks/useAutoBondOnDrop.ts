/**
 * useAutoBondOnDrop.ts
 * Hook that automatically creates bonds when atoms are added and autoBond is enabled.
 * Uses kernel bond validation logic.
 * 
 * This hook monitors the molecule for new atoms and automatically creates bonds
 * when autoBond is enabled, using the kernel's valence and distance validation.
 */

import { useEffect, useRef } from 'react';
import { shouldFormBond } from '../kernel/valence';
import { createBondSafe } from '../kernel/bonds';
import { useMoleculeStore } from '../store/moleculeStore';
import { pushState } from '../store/historyStore';

export function useAutoBondOnDrop() {
  const currentMolecule = useMoleculeStore((state) => state.currentMolecule);
  const autoBond = useMoleculeStore((state) => state.autoBond);
  const previousAtomCount = useRef(0);
  const previousMoleculeRef = useRef<typeof currentMolecule>(null);

  useEffect(() => {
    if (!autoBond || !currentMolecule) {
      previousAtomCount.current = currentMolecule?.atoms.size || 0;
      previousMoleculeRef.current = currentMolecule;
      return;
    }

    const currentAtomCount = currentMolecule.atoms.size;
    
    // Check if a new atom was added (atom count increased)
    if (currentAtomCount > previousAtomCount.current && previousAtomCount.current > 0) {
      const atoms = Array.from(currentMolecule.atoms.values());
      const latestAtom = atoms[atoms.length - 1];
      
      // Get all bonds as array for valence checking
      const bondsArray = Array.from(currentMolecule.bonds.values());
      
      // Try to bond with the closest atom that meets criteria
      let bonded = false;
      for (const atom of atoms.slice(0, -1)) {
        // Check if bond already exists
        const existingBond = bondsArray.find(
          (bond) => 
            (bond.a1 === atom.id && bond.a2 === latestAtom.id) ||
            (bond.a1 === latestAtom.id && bond.a2 === atom.id)
        );
        
        if (!existingBond && shouldFormBond(atom, latestAtom, bondsArray)) {
          // Create bond using kernel function
          const bondId = createBondSafe(
            currentMolecule,
            atom.id,
            latestAtom.id
          );
          
          if (bondId) {
            // Update store with new molecule state
            const cloned = currentMolecule.clone();
            useMoleculeStore.getState().setMolecule(cloned);
            pushState();
            bonded = true;
            break; // Only bond to one atom at a time
          }
        }
      }
      
      // Update refs after processing
      if (bonded) {
        previousAtomCount.current = currentMolecule.atoms.size;
        previousMoleculeRef.current = currentMolecule;
      }
    }
    
    // Update refs
    previousAtomCount.current = currentAtomCount;
    previousMoleculeRef.current = currentMolecule;
  }, [currentMolecule, autoBond]);
}

