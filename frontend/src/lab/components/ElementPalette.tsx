/**
 * ElementPalette - Element selection palette
 */

import React from 'react';
import { useLab } from '../hooks/useLab';

const COMMON_ELEMENTS = [
  { symbol: 'C', name: 'Carbon', color: '#909090' },
  { symbol: 'H', name: 'Hydrogen', color: '#ffffff' },
  { symbol: 'O', name: 'Oxygen', color: '#ff0d0d' },
  { symbol: 'N', name: 'Nitrogen', color: '#3050f8' },
  { symbol: 'F', name: 'Fluorine', color: '#90e050' },
  { symbol: 'Cl', name: 'Chlorine', color: '#1ff01f' },
  { symbol: 'Br', name: 'Bromine', color: '#a62929' },
  { symbol: 'I', name: 'Iodine', color: '#940094' },
  { symbol: 'S', name: 'Sulfur', color: '#ffff30' },
  { symbol: 'P', name: 'Phosphorus', color: '#ff8000' },
  { symbol: 'B', name: 'Boron', color: '#ffb5b5' },
  { symbol: 'Si', name: 'Silicon', color: '#f0c8a0' },
];

export default function ElementPalette() {
  const { currentElement, setCurrentElement } = useLab();

  return (
    <div className="p-3 border-b border-neutral-700">
      <h2 className="text-sm font-semibold mb-2 opacity-70">Elements</h2>
      <div className="grid grid-cols-3 gap-2">
        {COMMON_ELEMENTS.map((element) => (
          <button
            key={element.symbol}
            onClick={() => setCurrentElement(element.symbol)}
            className={`px-2 py-2 rounded text-xs font-semibold transition-all ${
              currentElement === element.symbol
                ? 'bg-blue-600 text-white ring-2 ring-blue-400'
                : 'bg-neutral-700 text-neutral-300 hover:bg-neutral-600'
            }`}
            style={{
              backgroundColor:
                currentElement === element.symbol ? undefined : element.color + '20',
              borderColor: element.color,
            }}
            title={element.name}
          >
            <div className="flex flex-col items-center">
              <div
                className="w-4 h-4 rounded-full mb-1"
                style={{ backgroundColor: element.color }}
              />
              <span>{element.symbol}</span>
            </div>
          </button>
        ))}
      </div>
      <div className="mt-2 text-xs opacity-60">
        Selected: <span className="font-mono">{currentElement}</span>
      </div>
    </div>
  );
}

