import React from 'react';
import { motion } from 'framer-motion';
import Card from '../components/ui/Card';

interface ModelCard {
	id: string;
	name: string;
	description: string;
	version: string;
	endpoint: string;
	metrics?: {
		accuracy?: string;
		latency?: string;
	};
}

const models: ModelCard[] = [
	{
		id: 'coregen-1',
		name: 'CoreGen-1',
		description: 'Generates 3D molecule structures from text descriptions or SMILES notation. Powers the molecule generation pipeline.',
		version: 'v2.1.0',
		endpoint: '/api/molecule/generate',
		metrics: {
			accuracy: '94.2%',
			latency: '120ms',
		},
	},
	{
		id: 'propnet-x',
		name: 'PropNet-X',
		description: 'QSAR property predictor. Estimates solubility, pKa, toxicity, drug-likeness, and other molecular properties.',
		version: 'v1.8.3',
		endpoint: '/api/analysis/qsar',
		metrics: {
			accuracy: '91.5%',
			latency: '85ms',
		},
	},
	{
		id: 'reactflow-r1',
		name: 'ReactFlow-R1',
		description: 'Multi-step synthesis route planner. Suggests optimal reaction pathways for target molecules.',
		version: 'v1.2.0',
		endpoint: '/api/reaction/plan',
		metrics: {
			accuracy: '87.3%',
			latency: '450ms',
		},
	},
	{
		id: 'chemgpt-s',
		name: 'ChemGPT-S',
		description: 'Reaction simulator and outcome predictor. Models chemical reactions and predicts products with high fidelity.',
		version: 'v1.5.2',
		endpoint: '/api/reaction/simulate',
		metrics: {
			accuracy: '89.7%',
			latency: '200ms',
		},
	},
];

export default function Models() {
	return (
		<motion.div
			initial={{ opacity: 0 }}
			animate={{ opacity: 1 }}
			transition={{ duration: 0.3 }}
			className="space-y-8"
		>
			{/* Header */}
			<div className="space-y-2">
				<h1 className="text-4xl font-bold text-black">AI Models</h1>
				<p className="text-darkGrey text-lg">
					Explore our machine learning models for molecular design and analysis
				</p>
			</div>

			{/* Model Grid */}
			<div className="grid grid-cols-1 md:grid-cols-2 gap-6">
				{models.map((model) => (
					<Card key={model.id} className="p-6 hover:shadow-neon-hover transition-shadow">
						<div className="space-y-4">
							<div className="flex items-start justify-between">
								<div>
									<h3 className="text-2xl font-semibold text-black">{model.name}</h3>
									<p className="text-sm text-midGrey mt-1">Version {model.version}</p>
								</div>
							</div>
							
							<p className="text-darkGrey leading-relaxed">{model.description}</p>
							
							{model.metrics && (
								<div className="flex gap-4 pt-2 border-t border-lightGrey">
									{model.metrics.accuracy && (
										<div>
											<span className="text-xs text-midGrey">Accuracy</span>
											<div className="text-lg font-semibold text-black">{model.metrics.accuracy}</div>
										</div>
									)}
									{model.metrics.latency && (
										<div>
											<span className="text-xs text-midGrey">Latency</span>
											<div className="text-lg font-semibold text-black">{model.metrics.latency}</div>
										</div>
									)}
								</div>
							)}
							
							<div className="pt-2">
								<code className="text-xs bg-offwhite px-2 py-1 rounded border border-lightGrey text-darkGrey">
									{model.endpoint}
								</code>
							</div>
							
							<div className="pt-2">
								<button className="btn-secondary px-4 py-2 text-sm">
									Try it out
								</button>
							</div>
						</div>
					</Card>
				))}
			</div>

			{/* Additional Info */}
			<Card className="p-6 bg-offwhite">
				<h3 className="text-xl font-semibold text-black mb-3">Model Access</h3>
				<p className="text-darkGrey mb-4">
					All models are accessible via REST API endpoints. Check the{' '}
					<a href="/docs" className="text-black underline hover:text-darkGrey">
						Documentation
					</a>{' '}
					for detailed API reference, authentication, and usage examples.
				</p>
				<div className="flex gap-3">
					<button className="btn-primary px-4 py-2 text-sm">
						View API Docs
					</button>
					<button className="btn-secondary px-4 py-2 text-sm">
						Get API Key
					</button>
				</div>
			</Card>
		</motion.div>
	);
}

