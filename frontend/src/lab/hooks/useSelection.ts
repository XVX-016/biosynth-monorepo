/**
 * useSelection - Hook for accessing selection state
 */

import { useLabStore } from '../state/LabStore';

export function useSelection() {
  const selectedAtomId = useLabStore((s) => s.selectedAtomId);
  const selectedBondId = useLabStore((s) => s.selectedBondId);
  const selectAtom = useLabStore((s) => s.selectAtom);
  const selectBond = useLabStore((s) => s.selectBond);
  
  return {
    selectedAtomId,
    selectedBondId,
    selectAtom,
    selectBond,
    clear: () => {
      selectAtom(null);
      selectBond(null);
    },
  };
}

