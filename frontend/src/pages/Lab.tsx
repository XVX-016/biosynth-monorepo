import React, { useState } from 'react';
import MoleculeViewer from '../components/MoleculeViewer';
import ToolPanel from '../components/ToolPanel';
import PropertiesPanel from '../components/PropertiesPanel';
import { TemplatePanel } from '../components/TemplatePanel';
import Button from '../components/ui/Button';
import { useMoleculeStore } from '../store/moleculeStore';
import { saveMolecule } from '../lib/api';
import { moleculeToJSON, getCanvasThumbnail } from '../lib/engineAdapter';
import { MoleculeSerializer } from '@biosynth/engine';

export default function Lab() {
	const molecule = useMoleculeStore((s) => s.currentMolecule);
	const backendPredictions = useMoleculeStore((s) => s.backendPredictions);
	const [saving, setSaving] = useState(false);

	const handleSave = async () => {
		if (!molecule) return;
		const name = prompt('Enter molecule name:');
		if (!name) return;
		setSaving(true);
		try {
			const jsonGraph = moleculeToJSON(molecule);
			const smiles = MoleculeSerializer.toSMILES(molecule);
			const thumbnail = getCanvasThumbnail();
			const properties = backendPredictions ? JSON.stringify(backendPredictions) : undefined;
			await saveMolecule({
				name,
				smiles: smiles || undefined,
				json_graph: jsonGraph,
				properties,
				thumbnail_b64: thumbnail,
			});
			alert(`Molecule "${name}" saved.`);
		} catch (e) {
			alert('Failed to save molecule');
		} finally {
			setSaving(false);
		}
	};

	return (
		<div className="lab-layout">
			<TemplatePanel />
			<div className="grid grid-cols-12 gap-4 flex-1">
				<div className="col-span-12 lg:col-span-3">
					<div className="bg-white rounded-xl shadow-neon border border-lightGrey p-4">
						<ToolPanel />
					</div>
				</div>
				<div className="col-span-12 lg:col-span-6">
					<div className="relative bg-white rounded-xl shadow-neon border border-lightGrey p-2 h-[70vh] lg:h-[78vh] bg-offwhite">
						<MoleculeViewer />
						<div className="absolute bottom-3 right-3">
							<Button onClick={handleSave} disabled={!molecule || saving}>
								{saving ? 'Saving...' : 'Save'}
							</Button>
						</div>
					</div>
				</div>
				<div className="col-span-12 lg:col-span-3">
					<div className="bg-white rounded-xl shadow-neon border border-lightGrey p-4">
						<PropertiesPanel />
					</div>
				</div>
			</div>
		</div>
	);
}


