/**
 * ValidationPanel - Enhanced validation display with auto-fix suggestions
 */

import React from 'react';
import { useLab } from '../hooks/useLab';
import { useLabStore } from '../state/LabStore';

export default function ValidationPanel() {
  const { validationResult, moleculeEngine, dispatch } = useLab();

  const handleAutoFix = (errorType: string) => {
    if (!validationResult) return;

    switch (errorType) {
      case 'duplicate_bond':
        // Remove duplicate bonds
        const seen = new Set<string>();
        const bonds = Array.from(moleculeEngine.getBonds());
        bonds.forEach((bond) => {
          const key = [bond.atoms[0], bond.atoms[1]].sort().join('-');
          if (seen.has(key)) {
            dispatch('deleteBond', {
              bondId: bond.id,
              bond,
            });
          }
          seen.add(key);
        });
        break;

      case 'charge':
        // Reset unrealistic charges
        moleculeEngine.getAtoms().forEach((atom) => {
          if (Math.abs(atom.charge) > 3) {
            // Would need an updateAtom action
            console.log('Would reset charge for atom', atom.id);
          }
        });
        break;

      default:
        console.log('Auto-fix not implemented for', errorType);
    }
  };

  const handleSanitize = () => {
    const { validationEngine } = useLabStore.getState();
    const state = moleculeEngine.getState();
    const sanitized = validationEngine.sanitize(state);
    // Reload sanitized state
    moleculeEngine.loadState(sanitized);
    // Re-validate and update store
    useLabStore.getState().validate();
  };

  if (!validationResult) {
    return (
      <div className="p-3">
        <h2 className="text-sm font-semibold mb-2 opacity-70">Validation</h2>
        <div className="text-xs opacity-60">No validation performed</div>
      </div>
    );
  }

  const errorsByType = validationResult.errors.reduce(
    (acc, error) => {
      if (!acc[error.type]) acc[error.type] = [];
      acc[error.type].push(error);
      return acc;
    },
    {} as Record<string, typeof validationResult.errors>
  );

  const canAutoFix = Object.keys(errorsByType).some(
    (type) => type === 'duplicate_bond' || type === 'charge'
  );

  return (
    <div className="p-3 border-t border-neutral-700">
      <div className="flex items-center justify-between mb-2">
        <h2 className="text-sm font-semibold opacity-70">Validation</h2>
        {canAutoFix && (
          <button
            onClick={handleSanitize}
            className="text-xs px-2 py-1 bg-blue-600 text-white rounded hover:bg-blue-700"
          >
            Auto-Fix
          </button>
        )}
      </div>

      {validationResult.valid && validationResult.errors.length === 0 ? (
        <div className="text-green-400 text-xs">✓ Molecule is valid</div>
      ) : (
        <div className="space-y-2">
          {/* Errors by type */}
          {Object.entries(errorsByType).map(([type, errors]) => (
            <div key={type} className="bg-neutral-700 p-2 rounded">
              <div className="flex items-center justify-between mb-1">
                <div className="text-red-400 text-xs font-semibold">
                  {type} ({errors.length})
                </div>
                {(type === 'duplicate_bond' || type === 'charge') && (
                  <button
                    onClick={() => handleAutoFix(type)}
                    className="text-xs px-1 py-0.5 bg-blue-600 text-white rounded hover:bg-blue-700"
                  >
                    Fix
                  </button>
                )}
              </div>
              <div className="space-y-1">
                {errors.slice(0, 3).map((error, i) => (
                  <div key={i} className="text-xs text-red-300">
                    • {error.message}
                    {error.atomId && (
                      <span className="ml-1 opacity-60">(atom: {error.atomId.slice(0, 8)})</span>
                    )}
                    {error.bondId && (
                      <span className="ml-1 opacity-60">(bond: {error.bondId.slice(0, 8)})</span>
                    )}
                  </div>
                ))}
                {errors.length > 3 && (
                  <div className="text-xs opacity-60">+ {errors.length - 3} more...</div>
                )}
              </div>
            </div>
          ))}

          {/* Warnings */}
          {validationResult.warnings.length > 0 && (
            <div className="bg-neutral-700 p-2 rounded">
              <div className="text-yellow-400 text-xs font-semibold mb-1">
                Warnings ({validationResult.warnings.length})
              </div>
              <div className="space-y-1">
                {validationResult.warnings.map((warning, i) => (
                  <div key={i} className="text-xs text-yellow-300">
                    • {warning}
                  </div>
                ))}
              </div>
            </div>
          )}
        </div>
      )}
    </div>
  );
}

