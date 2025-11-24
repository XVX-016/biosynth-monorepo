import React from 'react'
import { Html } from '@react-three/drei'
import type { Atom } from '../../types/molecule'

interface AtomTooltipProps {
  atom: Atom
  position: [number, number, number]
}

export default function AtomTooltip({ atom, position }: AtomTooltipProps) {
  return (
    <Html position={position} center>
      <div className="px-2 py-1 rounded-md text-xs bg-white/90 backdrop-blur-sm border border-gray-300 shadow-sm">
        <div className="font-semibold">{atom.element}</div>
        <div className="text-gray-500 text-[10px]">ID: {atom.id.slice(0, 8)}</div>
      </div>
    </Html>
  )
}

