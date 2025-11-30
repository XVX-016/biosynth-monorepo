/**
 * useLab - Main hook for accessing Lab store
 */

import { useLabStore } from '../state/LabStore';

export function useLab() {
  return useLabStore();
}

