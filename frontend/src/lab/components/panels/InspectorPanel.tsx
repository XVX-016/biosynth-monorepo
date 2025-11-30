/**
 * InspectorPanel - Panel for inspecting selected atoms/bonds
 */

import React from 'react';
import { useSelection } from '../../hooks/useSelection';
import { useLab } from '../../hooks/useLab';

export default function InspectorPanel() {
  const { selectedAtomId, selectedBondId } = useSelection();
  const { moleculeEngine } = useLab();

  const selectedAtom = selectedAtomId
    ? moleculeEngine.getAtom(selectedAtomId)
    : null;
  const selectedBond = selectedBondId
    ? moleculeEngine.getBond(selectedBondId)
    : null;

  return (
    <div className="p-3">
      <h2 className="text-sm font-semibold mb-2 opacity-70">Inspector</h2>
      
      {selectedAtom && (
        <div className="bg-neutral-700 p-2 rounded text-sm">
          <div className="font-semibold mb-1">Atom</div>
          <div>Element: {selectedAtom.element}</div>
          <div>Position: ({selectedAtom.x.toFixed(2)}, {selectedAtom.y.toFixed(2)})</div>
          <div>Charge: {selectedAtom.charge}</div>
        </div>
      )}
      
      {selectedBond && (
        <div className="bg-neutral-700 p-2 rounded text-sm mt-2">
          <div className="font-semibold mb-1">Bond</div>
          <div>Order: {selectedBond.order}</div>
          <div>Atoms: {selectedBond.atoms.join(' - ')}</div>
        </div>
      )}
      
      {!selectedAtom && !selectedBond && (
        <div className="text-xs opacity-60">No selection</div>
      )}
    </div>
  );
}

