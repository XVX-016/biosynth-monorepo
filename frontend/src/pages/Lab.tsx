import React from 'react';
import MoleculeViewer from '../components/MoleculeViewer';
import ToolPanel from '../components/ToolPanel';

export default function Lab() {
	return (
		<div className="grid grid-cols-12 gap-4">
			<div className="col-span-12 lg:col-span-3">
				<div className="bg-panel rounded-xl shadow-soft border border-aluminum-DEFAULT p-4">
					<ToolPanel />
				</div>
			</div>
			<div className="col-span-12 lg:col-span-6">
				<div className="bg-panel rounded-xl shadow-soft border border-aluminum-DEFAULT p-2 h-[70vh] lg:h-[78vh]">
					<MoleculeViewer />
				</div>
			</div>
			<div className="col-span-12 lg:col-span-3">
				<div className="bg-panel rounded-xl shadow-soft border border-aluminum-DEFAULT p-4">
					<h3 className="text-lg font-semibold mb-2">Properties</h3>
					<p className="text-text-secondary">Select an atom or bond to see details.</p>
				</div>
			</div>
		</div>
	);
}


