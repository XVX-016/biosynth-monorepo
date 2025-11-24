import React from 'react'
import { useLabStore } from '../../store/labStore'
import elementPalette from '../../data/elementPalette'

const ELEMENTS = ['C', 'H', 'O', 'N', 'S', 'P', 'F', 'Cl', 'Br', 'I'] as const

export default function ElementPalette() {
  const current = useLabStore(s => s.currentElement)
  const set = useLabStore(s => s.setCurrentElement)
  
  return (
    <div className="grid grid-cols-3 gap-2 w-full">
      {ELEMENTS.map(sym => {
        const pal = elementPalette[sym] || elementPalette.C
        const baseColor = `#${pal.base.toString(16).padStart(6, '0')}`
        const accentColor = `#${pal.accent.toString(16).padStart(6, '0')}`
        const isSelected = current === sym
        
        return (
          <button
            key={sym}
            onClick={() => set(sym)}
            className={`
              w-full aspect-square rounded-md border transition-all
              flex items-center justify-center text-xs font-bold
              ${isSelected 
                ? 'border-[#4676ff] ring-2 ring-[#4676ff]/20 shadow-sm' 
                : 'border-gray-300 hover:border-gray-400'}
            `}
            style={{
              background: baseColor,
              color: accentColor,
            }}
            title={sym}
          >
            {sym}
          </button>
        )
      })}
    </div>
  )
}

