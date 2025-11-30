import { memo } from 'react'
import type { Element } from '@biosynth/engine'
import { useMoleculeStore } from '../store/moleculeStore'

const ELEMENTS: Element[] = ['C', 'H', 'O', 'N', 'S', 'P', 'F', 'Cl']

/**
 * Minimal element palette tied to the molecule store.
 * Selecting an element auto-enables the add-atom tool for quicker placement.
 */
function LabElementPalette() {
  const activeElement = useMoleculeStore((s) => s.atomToAdd)
  const setAtomToAdd = useMoleculeStore((s) => s.setAtomToAdd)
  const setTool = useMoleculeStore((s) => s.setTool)

  const handleSelect = (element: Element) => {
    setAtomToAdd(element)
    setTool('add-atom')
  }

  return (
    <div className="grid grid-cols-4 gap-2">
      {ELEMENTS.map((element) => {
        const isActive = activeElement === element
        return (
          <button
            key={element}
            onClick={() => handleSelect(element)}
            className={`rounded-lg border px-0 py-3 text-sm font-semibold transition-all ${
              isActive
                ? 'bg-black text-white border-black shadow-neon-sm'
                : 'bg-white text-darkGrey border-lightGrey hover:border-darkGrey/40 hover:text-black'
            }`}
            aria-pressed={isActive}
          >
            {element}
          </button>
        )
      })}
    </div>
  )
}

export default memo(LabElementPalette)


