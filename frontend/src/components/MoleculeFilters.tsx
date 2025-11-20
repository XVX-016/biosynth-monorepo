/**
 * MoleculeFilters - Search and filter component for molecule library
 * 
 * Supports:
 * - Text search (name, formula, SMILES)
 * - Element filter
 * - Debounced input
 */
import React, { useState, useEffect } from 'react';

export interface MoleculeFilters {
  query: string;
  element: string;
}

interface MoleculeFiltersProps {
  value?: MoleculeFilters;
  onChange: (filters: MoleculeFilters) => void;
  placeholder?: string;
  debounceMs?: number;
}

export default function MoleculeFilters({
  value = { query: '', element: '' },
  onChange,
  placeholder = 'Search by name, formula or SMILES...',
  debounceMs = 350,
}: MoleculeFiltersProps) {
  const [query, setQuery] = useState(value.query);
  const [element, setElement] = useState(value.element);

  // Sync with external value changes
  useEffect(() => {
    setQuery(value.query);
    setElement(value.element);
  }, [value]);

  // Debounce onChange callback
  useEffect(() => {
    const timer = setTimeout(() => {
      onChange({ query: query.trim(), element: element.trim() });
    }, debounceMs);

    return () => clearTimeout(timer);
  }, [query, element, onChange, debounceMs]);

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    onChange({ query: query.trim(), element: element.trim() });
  };

  return (
    <form onSubmit={handleSubmit} className="flex gap-3 items-center">
      <input
        type="text"
        value={query}
        onChange={(e) => setQuery(e.target.value)}
        placeholder={placeholder}
        className="flex-1 rounded-full border border-gray-200 px-4 py-3 shadow-sm focus:outline-none focus:ring-2 focus:ring-darkGrey/20 focus:border-darkGrey placeholder:text-midGrey"
      />
      <input
        type="text"
        value={element}
        onChange={(e) => setElement(e.target.value.toUpperCase())}
        placeholder="Element (C, O, N)"
        maxLength={2}
        className="w-32 rounded-full border border-gray-200 px-4 py-3 shadow-sm focus:outline-none focus:ring-2 focus:ring-darkGrey/20 focus:border-darkGrey placeholder:text-midGrey text-center font-mono"
      />
      <button
        type="submit"
        className="px-4 py-2 rounded-full bg-black text-white hover:bg-darkGrey transition-colors whitespace-nowrap"
      >
        Search
      </button>
      {(query || element) && (
        <button
          type="button"
          onClick={() => {
            setQuery('');
            setElement('');
            onChange({ query: '', element: '' });
          }}
          className="px-3 py-2 rounded-full bg-gray-200 text-gray-700 hover:bg-gray-300 transition-colors text-sm"
        >
          Clear
        </button>
      )}
    </form>
  );
}

