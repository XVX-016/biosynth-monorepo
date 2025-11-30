/**
 * BottomDock - Bottom panel for validation, logs, etc.
 */

import React from 'react';
import { useLab } from '../hooks/useLab';

export default function BottomDock() {
  const { validationResult, currentMolecule } = useLab();

  return (
    <div className="h-32 border-t border-neutral-700 bg-neutral-800 p-3 overflow-y-auto">
      <div className="text-xs font-semibold mb-2 opacity-70">Validation</div>
      
      {validationResult && (
        <div className="space-y-1">
          {validationResult.errors.length > 0 && (
            <div className="text-red-400 text-xs">
              {validationResult.errors.length} error(s):
              {validationResult.errors.map((err, i) => (
                <div key={i} className="ml-2">• {err.message}</div>
              ))}
            </div>
          )}
          {validationResult.warnings.length > 0 && (
            <div className="text-yellow-400 text-xs">
              {validationResult.warnings.length} warning(s):
              {validationResult.warnings.map((warn, i) => (
                <div key={i} className="ml-2">• {warn}</div>
              ))}
            </div>
          )}
          {validationResult.valid && validationResult.errors.length === 0 && (
            <div className="text-green-400 text-xs">✓ Molecule is valid</div>
          )}
        </div>
      )}
      
      {!validationResult && (
        <div className="text-xs opacity-60">No validation performed</div>
      )}
    </div>
  );
}

