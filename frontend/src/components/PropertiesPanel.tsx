import React, { useMemo } from 'react';
import { useMoleculeStore } from '../store/moleculeStore';

export default function PropertiesPanel() {
	const molecule = useMoleculeStore((s) => s.currentMolecule);
	const selectedAtomId = useMoleculeStore((s) => s.selectedAtomId);
	const selectedBondId = useMoleculeStore((s) => s.selectedBondId);

	const content = useMemo(() => {
		if (!molecule) {
			return <div className="text-chrome">No molecule loaded.</div>;
		}
		if (selectedAtomId) {
			const atom = molecule.atoms.get(selectedAtomId);
			if (!atom) return null;
			return (
				<div>
					<div className="text-sm text-chrome mb-1">Selected Atom</div>
					<div className="text-xl font-semibold text-ivory">{atom.element}</div>
					<div className="text-sm text-chrome mt-2">
						Position: ({atom.position[0].toFixed(2)}, {atom.position[1].toFixed(2)},{' '}
						{atom.position[2].toFixed(2)})
					</div>
					<div className="text-sm text-chrome/70 mt-1">ID: {atom.id}</div>
				</div>
			);
		}
		if (selectedBondId) {
			const bond = molecule.bonds.get(selectedBondId);
			if (!bond) return null;
			const a1 = molecule.atoms.get(bond.a1);
			const a2 = molecule.atoms.get(bond.a2);
			return (
				<div>
					<div className="text-sm text-chrome mb-1">Selected Bond</div>
					<div className="text-xl font-semibold text-ivory">
						{bond.order === 1 ? 'Single' : bond.order === 2 ? 'Double' : 'Triple'}
					</div>
					<div className="text-sm text-chrome mt-2">
						Between: {a1?.element} - {a2?.element}
					</div>
					<div className="text-sm text-chrome/70 mt-1">ID: {bond.id}</div>
				</div>
			);
		}
		return <div className="text-chrome">Select an atom or bond to view details.</div>;
	}, [molecule, selectedAtomId, selectedBondId]);

	return (
		<div className="space-y-4">
			<h3 className="text-lg font-semibold text-ivory">Properties</h3>
			{content}
		</div>
	);
}


