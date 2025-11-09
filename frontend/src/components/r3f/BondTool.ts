import { useEffect, useRef } from 'react'
import { useMoleculeStore } from '../../store/moleculeStore'
import { updateAtomPosition, addBond } from '../../lib/engineAdapter'

/**
 * BondTool - Handles bond creation between selected atoms
 * Creates a bond when two atoms are selected within 2 seconds
 */
export function useBondTool() {
  const selectedAtomId = useMoleculeStore((state) => state.selectedAtomId)
  const currentMolecule = useMoleculeStore((state) => state.currentMolecule)
  const firstSelectedRef = useRef<string | null>(null)
  const timeoutRef = useRef<NodeJS.Timeout | null>(null)

  useEffect(() => {
    // Reset timeout when selection changes
    if (timeoutRef.current) {
      clearTimeout(timeoutRef.current)
      timeoutRef.current = null
    }

    if (!selectedAtomId || !currentMolecule) {
      firstSelectedRef.current = null
      return
    }

    // If no first atom selected, store it
    if (!firstSelectedRef.current) {
      firstSelectedRef.current = selectedAtomId
      return
    }

    // If same atom selected, ignore
    if (firstSelectedRef.current === selectedAtomId) {
      return
    }

    // Two different atoms selected - create bond
    const atom1 = firstSelectedRef.current
    const atom2 = selectedAtomId

    // Check if bond already exists
    const existingBond = Array.from(currentMolecule.bonds.values()).find(
      (b) =>
        (b.a1 === atom1 && b.a2 === atom2) ||
        (b.a1 === atom2 && b.a2 === atom1)
    )

    if (!existingBond) {
      // Create bond
      addBond(atom1, atom2)
    }

    // Reset after bond creation
    firstSelectedRef.current = null
  }, [selectedAtomId, currentMolecule])

  // Handle Escape key to cancel bond creation
  useEffect(() => {
    const handleEscape = (e: KeyboardEvent) => {
      if (e.key === 'Escape') {
        firstSelectedRef.current = null
        if (timeoutRef.current) {
          clearTimeout(timeoutRef.current)
          timeoutRef.current = null
        }
      }
    }

    window.addEventListener('keydown', handleEscape)
    return () => {
      window.removeEventListener('keydown', handleEscape)
    }
  }, [])

  return {
    firstSelectedAtomId: firstSelectedRef.current,
    isWaitingForSecondAtom: firstSelectedRef.current !== null,
  }
}

