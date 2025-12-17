import React, { useState } from 'react';
import { motion } from 'framer-motion';
import Card from '../components/ui/Card';

interface DocSection {
	id: string;
	title: string;
	content: React.ReactNode;
}

const sections: DocSection[] = [
	{
		id: 'getting-started',
		title: 'Getting Started',
		content: (
			<div className="space-y-4">
				<p className="text-darkGrey">
					Welcome to MolForge documentation. This guide will help you get started with our molecular design platform.
				</p>
				<h3 className="text-xl font-semibold text-black mt-6">Quick Start</h3>
				<ol className="list-decimal list-inside space-y-2 text-darkGrey ml-4">
					<li>Create an account or sign in</li>
					<li>Navigate to the Lab to start designing molecules</li>
					<li>Use the Library to browse and manage your saved molecules</li>

				</ol>
			</div>
		),
	},
	{
		id: 'using-lab',
		title: 'Using the Lab',
		content: (
			<div className="space-y-4">
				<p className="text-darkGrey">
					The Lab is your interactive molecular editor. Build molecules atom by atom, adjust bonds, and visualize structures in 3D.
				</p>
				<h3 className="text-xl font-semibold text-black mt-6">Features</h3>
				<ul className="list-disc list-inside space-y-2 text-darkGrey ml-4">
					<li>Drag-and-drop atom placement</li>
					<li>Automatic bond formation</li>
					<li>3D visualization with rotation and zoom</li>
					<li>Export to MOL, PDB, SMILES, PNG, and GLTF formats</li>
					<li>Real-time property predictions</li>
				</ul>
			</div>
		),
	},
	{
		id: 'using-library',
		title: 'Using the Library',
		content: (
			<div className="space-y-4">
				<p className="text-darkGrey">
					The Library stores all your saved molecules. Search, filter, and organize your molecular collection.
				</p>
				<h3 className="text-xl font-semibold text-black mt-6">Features</h3>
				<ul className="list-disc list-inside space-y-2 text-darkGrey ml-4">
					<li>Search by name or SMILES notation</li>
					<li>Filter by molecule type</li>
					<li>3D preview thumbnails</li>
					<li>Load molecules directly into the Lab</li>
					<li>Export and share molecules</li>
				</ul>
			</div>
		),
	},
	{
		id: 'api-reference',
		title: 'API Reference',
		content: (
			<div className="space-y-6">
				<p className="text-darkGrey">
					MolForge provides a RESTful API for programmatic access to all features.
				</p>

				<div className="space-y-4">
					<div>
						<h4 className="text-lg font-semibold text-black mb-2">POST /api/molecule/generate</h4>
						<p className="text-darkGrey mb-2">Generate a 3D molecule structure from SMILES or text description.</p>
						<pre className="bg-offwhite p-3 rounded border border-lightGrey text-sm text-darkGrey overflow-x-auto">
							{`{
  "input": "CCO",
  "format": "smiles"
}`}
						</pre>
					</div>

					<div>
						<h4 className="text-lg font-semibold text-black mb-2">POST /api/analysis/qsar</h4>
						<p className="text-darkGrey mb-2">Predict molecular properties using QSAR models.</p>
						<pre className="bg-offwhite p-3 rounded border border-lightGrey text-sm text-darkGrey overflow-x-auto">
							{`{
  "smiles": "CCO",
  "properties": ["solubility", "toxicity"]
}`}
						</pre>
					</div>

					<div>
						<h4 className="text-lg font-semibold text-black mb-2">POST /api/reaction/simulate</h4>
						<p className="text-darkGrey mb-2">Simulate a chemical reaction and predict products.</p>
						<pre className="bg-offwhite p-3 rounded border border-lightGrey text-sm text-darkGrey overflow-x-auto">
							{`{
  "reactants": ["CCO", "O"],
  "conditions": {...}
}`}
						</pre>
					</div>
				</div>
			</div>
		),
	},
	{
		id: 'data-formats',
		title: 'Data Formats',
		content: (
			<div className="space-y-4">
				<p className="text-darkGrey">
					MolForge supports multiple molecular data formats for import and export.
				</p>

				<div className="space-y-3">
					<div>
						<h4 className="font-semibold text-black">SMILES</h4>
						<p className="text-darkGrey text-sm">Simplified Molecular Input Line Entry System - text-based representation.</p>
						<code className="text-xs bg-offwhite px-2 py-1 rounded border border-lightGrey text-darkGrey">CCO</code>
					</div>

					<div>
						<h4 className="font-semibold text-black">MOL</h4>
						<p className="text-darkGrey text-sm">MDL Molfile format - standard for molecular structure data.</p>
					</div>

					<div>
						<h4 className="font-semibold text-black">PDB</h4>
						<p className="text-darkGrey text-sm">Protein Data Bank format - commonly used for protein structures.</p>
					</div>

					<div>
						<h4 className="font-semibold text-black">GLTF</h4>
						<p className="text-darkGrey text-sm">3D model format for visualization and sharing.</p>
					</div>
				</div>
			</div>
		),
	},
	{
		id: 'viewer-embeds',
		title: 'Viewer Embeds',
		content: (
			<div className="space-y-4">
				<p className="text-darkGrey">
					Embed the MolForge 3D molecule viewer in your own applications.
				</p>
				<pre className="bg-offwhite p-3 rounded border border-lightGrey text-sm text-darkGrey overflow-x-auto">
					{`<iframe
  src="https://molforge.app/viewer?smiles=CCO"
  width="600"
  height="400"
  frameborder="0"
></iframe>`}
				</pre>
			</div>
		),
	},
];

export default function Docs() {
	const [activeSection, setActiveSection] = useState(sections[0].id);

	return (
		<motion.div
			initial={{ opacity: 0 }}
			animate={{ opacity: 1 }}
			transition={{ duration: 0.3 }}
			className="max-w-4xl mx-auto space-y-8"
		>
			{/* Header */}
			<motion.div
				className="space-y-2"
				initial={{ opacity: 0, y: -20 }}
				animate={{ opacity: 1, y: 0 }}
				transition={{ duration: 0.4, delay: 0.1 }}
			>
				<h1 className="text-4xl font-bold text-black">Documentation</h1>
				<p className="text-darkGrey text-lg">
					Complete guide to using MolForge
				</p>
			</motion.div>

			{/* Section Navigation */}
			<motion.div
				className="flex flex-wrap gap-2 border-b border-lightGrey pb-4"
				initial={{ opacity: 0, y: 10 }}
				animate={{ opacity: 1, y: 0 }}
				transition={{ duration: 0.4, delay: 0.2 }}
			>
				{sections.map((section) => (
					<button
						key={section.id}
						onClick={() => setActiveSection(section.id)}
						className={`px-4 py-2 rounded-lg font-medium transition-all duration-200 ${activeSection === section.id
								? 'bg-black text-white'
								: 'text-darkGrey hover:text-black hover:bg-offwhite'
							}`}
					>
						{section.title}
					</button>
				))}
			</motion.div>

			{/* Active Section Content */}
			<motion.div
				key={activeSection}
				initial={{ opacity: 0, y: 20 }}
				animate={{ opacity: 1, y: 0 }}
				transition={{ duration: 0.4, delay: 0.1 }}
			>
				<Card className="p-8">
					{sections.find((s) => s.id === activeSection)?.content}
				</Card>
			</motion.div>
		</motion.div>
	);
}

