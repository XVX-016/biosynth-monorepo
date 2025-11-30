import React, { useState, useEffect, useMemo } from 'react';
import { useSearchParams } from 'react-router-dom';
import { motion } from 'framer-motion';
import MoleculeViewer from '../components/MoleculeViewer';
import ToolPanel from '../components/ToolPanel';
import PropertiesPanel from '../components/PropertiesPanel';
import KABPanel from '../components/kab/KABPanel';
import { TemplatePanel } from '../components/TemplatePanel';
import Button from '../components/ui/Button';
import { useMoleculeStore } from '../store/moleculeStore';
import { saveUserMolecule, getUserMolecule } from '../lib/userMoleculeStore';
import { getPublicMolecule } from '../lib/publicMoleculeStore';
import { moleculeToJSON, getCanvasThumbnail, moleculeFromJSON } from '../lib/engineAdapter';
import { MoleculeSerializer, MoleculeGraph } from '@biosynth/engine';
import { supabase } from '../supabase';
import { convertSMILESToMolfile, generateThumbnailBase64 } from '../lib/api';
import LabElementPalette from '../components/LabElementPalette';
import { analyzeStructure } from '../utils/moleculeIssues';
import { runGnassPipeline, moleculeFingerprint } from '../services/gnassService';
import VisualizationPanel from '../components/VisualizationPanel';
import ExportPanel from '../components/ExportPanel';
import { decodeSharePayload } from '../utils/shareLink';

export default function Lab() {
	const [searchParams] = useSearchParams();
	const molecule = useMoleculeStore((s) => s.currentMolecule);
	const setMolecule = useMoleculeStore((s) => s.setMolecule);
	const backendPredictions = useMoleculeStore((s) => s.backendPredictions);
	const [saving, setSaving] = useState(false);
	const [userId, setUserId] = useState<string | null>(null);
	const [loadingMolecule, setLoadingMolecule] = useState(false);
	const [gnassStatus, setGnassStatus] = useState<'idle' | 'running' | 'ready' | 'error'>('idle');
	const [gnassSummary, setGnassSummary] = useState<Record<string, number> | null>(null);
	const [gnassError, setGnassError] = useState<string | null>(null);

	const moleculeStats = useMemo(() => {
		return {
			atoms: molecule?.atoms.size ?? 0,
			bonds: molecule?.bonds.size ?? 0,
			formula: calculateFormula(molecule),
		};
	}, [molecule]);

	const issueSummary = useMemo(() => analyzeStructure(molecule), [molecule]);
	const moleculeKey = useMemo(() => moleculeFingerprint(molecule), [molecule]);

	useEffect(() => {
		if (!supabase) return;

		// Check current session
		supabase.auth.getSession().then(({ data: { session } }) => {
			if (session?.user) {
				setUserId(session.user.id);
			}
		});

		// Listen for auth changes
		const {
			data: { subscription },
		} = supabase.auth.onAuthStateChange((_event, session) => {
			if (session?.user) {
				setUserId(session.user.id);
			} else {
				setUserId(null);
			}
		});

		return () => subscription.unsubscribe();
	}, []);

	// Load molecule from URL params
	useEffect(() => {
		const moleculeId = searchParams.get('id');
		const source = searchParams.get('source'); // 'public' or 'user'

		if (!moleculeId || !source) return;

		const loadMolecule = async () => {
			setLoadingMolecule(true);
			try {
				let molData;
				
				if (source === 'public') {
					molData = await getPublicMolecule(moleculeId);
				} else if (source === 'user' && userId) {
					molData = await getUserMolecule(userId, moleculeId);
				}

				if (!molData) {
					alert('Molecule not found');
					return;
				}

				// Try to load from metadata.json_graph first (preferred)
				if (molData.metadata?.json_graph) {
					const moleculeGraph = moleculeFromJSON(molData.metadata.json_graph);
					if (moleculeGraph) {
						setMolecule(moleculeGraph);
						setLoadingMolecule(false);
						return;
					}
				}

				// Fallback: Try to load from molfile or SMILES
				// Note: This would require converting molfile/SMILES to MoleculeGraph
				// For now, we'll just show an alert if json_graph is not available
				if (molData.molfile || molData.smiles) {
					alert('Molecule structure data is available but needs to be converted. Please use the template loader or create a new molecule.');
				} else {
					alert('Molecule data not available');
				}
			} catch (error) {
				console.error('Failed to load molecule:', error);
				alert('Failed to load molecule');
			} finally {
				setLoadingMolecule(false);
			}
		};

		// Wait for userId if loading user molecule
		if (source === 'user' && !userId) {
			// Wait a bit for auth to complete
			const timer = setTimeout(() => {
				if (userId) loadMolecule();
			}, 1000);
			return () => clearTimeout(timer);
		}

		loadMolecule();
	}, [searchParams, userId, setMolecule]);

	useEffect(() => {
		const shared = searchParams.get('share');
		if (!shared) return;
		try {
			const decoded = decodeSharePayload(shared);
			const sharedMolecule = moleculeFromJSON(decoded);
			if (sharedMolecule) {
				setMolecule(sharedMolecule);
				setLoadingMolecule(false);
			}
		} catch (error) {
			console.warn('Failed to import shared molecule', error);
		}
	}, [searchParams, setMolecule]);

	// Run GNASS/ML pipeline when structure changes
	useEffect(() => {
		if (!molecule || molecule.atoms.size === 0) {
			setGnassStatus('idle');
			setGnassSummary(null);
			setGnassError(null);
			return;
		}

		let cancelled = false;
		setGnassStatus('running');
		setGnassError(null);

		runGnassPipeline(molecule)
			.then((result) => {
				if (cancelled) return;
				setGnassSummary(result.predictions || null);
				setGnassStatus('ready');
			})
			.catch((error) => {
				if (cancelled) return;
				setGnassStatus('error');
				setGnassError(
					error instanceof Error ? error.message : 'GNASS pipeline failed. Retry after editing the molecule.',
				);
			});

		return () => {
			cancelled = true;
		};
		// eslint-disable-next-line react-hooks/exhaustive-deps
	}, [moleculeKey]);

	const handleSave = async () => {
		if (!molecule) return;
		
		if (!userId) {
			alert('Please sign in to save molecules. Visit /supabase-test to sign in.');
			return;
		}

		const name = prompt('Enter molecule name:');
		if (!name) return;
		
		setSaving(true);
		try {
			const jsonGraph = moleculeToJSON(molecule);
			const smiles = MoleculeSerializer.toSMILES(molecule);
			let thumbnail = getCanvasThumbnail();
			const properties = backendPredictions ? JSON.stringify(backendPredictions) : undefined;
			
			// Calculate formula if possible
			const formula = molecule ? calculateFormula(molecule) : undefined;
			
			// Generate molfile from SMILES if available
			let molfile: string | undefined = undefined;
			if (smiles && smiles.trim()) {
				try {
					const result = await convertSMILESToMolfile(smiles);
					molfile = result.molfile;
				} catch (e: any) {
					console.warn('Failed to generate molfile:', e);
					// Continue without molfile - it's optional
				}
			}
			
			// Generate thumbnail from backend if canvas thumbnail not available
			if (!thumbnail && smiles && smiles.trim()) {
				try {
					const thumbResult = await generateThumbnailBase64({ smiles });
					thumbnail = thumbResult.thumbnail_b64;
				} catch (e: any) {
					console.warn('Failed to generate thumbnail from backend:', e);
					// Continue without thumbnail - it's optional
				}
			}
			
			await saveUserMolecule(userId, {
				name,
				smiles: smiles || undefined,
				formula,
				molfile,
				metadata: {
					json_graph: jsonGraph,
					properties: properties ? JSON.parse(properties) : undefined,
				},
				thumbnail_b64: thumbnail || undefined,
			});
			alert(`Molecule "${name}" saved successfully!`);
		} catch (e: any) {
			console.error('Save error:', e);
			alert(`Failed to save molecule: ${e.message || 'Unknown error'}`);
		} finally {
			setSaving(false);
		}
	};

	// Helper function to calculate formula
	const calculateFormula = (mol: MoleculeGraph | null): string | undefined => {
		if (!mol) return undefined;
		const elementCounts: Record<string, number> = {};
		mol.atoms.forEach((atom) => {
			elementCounts[atom.element] = (elementCounts[atom.element] || 0) + 1;
		});
		return Object.entries(elementCounts)
			.map(([el, count]) => count > 1 ? `${el}${count}` : el)
			.join('');
	};

	return (
		<motion.div 
			className="lab-page min-h-screen bg-gradient-to-b from-offwhite via-white to-offwhite"
			initial={{ opacity: 0 }}
			animate={{ opacity: 1 }}
			transition={{ duration: 0.3 }}
		>
			{loadingMolecule && (
				<motion.div 
					className="fixed inset-0 bg-black/50 flex items-center justify-center z-50"
					initial={{ opacity: 0 }}
					animate={{ opacity: 1 }}
					transition={{ duration: 0.2 }}
				>
					<motion.div 
						className="bg-white rounded-lg p-6"
						initial={{ scale: 0.9, opacity: 0 }}
						animate={{ scale: 1, opacity: 1 }}
						transition={{ duration: 0.3 }}
					>
						<div className="text-center">Loading molecule...</div>
					</motion.div>
				</motion.div>
			)}
			<header className="lab-header">
				<div>
					<p className="text-sm uppercase tracking-[0.2em] text-midGrey">MolForge Lab</p>
					<h1 className="text-3xl font-semibold text-black mt-1">Interactive Molecule Workspace</h1>
				</div>
				<div className="lab-header__stats">
					<div>
						<p className="text-sm text-midGrey">Atoms</p>
						<p className="text-2xl font-semibold text-black">{moleculeStats.atoms}</p>
					</div>
					<div>
						<p className="text-sm text-midGrey">Bonds</p>
						<p className="text-2xl font-semibold text-black">{moleculeStats.bonds}</p>
					</div>
					<div>
						<p className="text-sm text-midGrey">Formula</p>
						<p className="text-xl font-semibold text-black">
							{moleculeStats.formula || (molecule ? 'Calculating…' : '—')}
						</p>
					</div>
				</div>
			</header>

			<div className="lab-grid">
				<section className="lab-sidebar space-y-6">
					<div className="card-base p-4 space-y-4">
						<div className="flex items-center justify-between">
							<h2 className="text-lg font-semibold text-black">Tools</h2>
							<span className="text-xs text-midGrey">Aligned controls</span>
						</div>
						<ToolPanel />
					</div>
					<div className="card-base p-4 space-y-3">
						<div className="flex items-center justify-between">
							<h3 className="text-base font-semibold text-black">Element Palette</h3>
							<span className="text-xs text-midGrey">Tap to place</span>
						</div>
						<LabElementPalette />
					</div>
					<div className="card-base p-4 space-y-3">
						<VisualizationPanel />
					</div>
					<div className="card-base p-4">
						<TemplatePanel />
					</div>
				</section>

				<section className="lab-workspace card-base">
					<div className="lab-workspace__canvas">
						<MoleculeViewer />
						<div className="lab-workspace__cta">
							<Button onClick={handleSave} disabled={!molecule || saving || !userId}>
								{saving ? 'Saving...' : userId ? 'Save to Library' : 'Sign in to Save'}
							</Button>
						</div>
					</div>
					<div className="lab-workspace__foot">
						<div>
							<p className="text-xs uppercase tracking-[0.3em] text-midGrey">Validation</p>
							<p className="text-sm text-black">{issueSummary.warnings[0]}</p>
						</div>
						<div className="text-right">
							<p className="text-xs uppercase tracking-[0.3em] text-midGrey">GNASS</p>
							<p className="text-sm text-black">
								{gnassStatus === 'running' && 'Training model…'}
								{gnassStatus === 'error' && 'Model error'}
								{gnassStatus === 'idle' && 'Awaiting molecule'}
								{gnassStatus === 'ready' && 'Predictions synced'}
							</p>
						</div>
					</div>
				</section>

				<section className="lab-inspector space-y-4">
					<div className="card-base p-4 space-y-4">
						<ExportPanel />
					</div>
					<div className="card-base p-4 space-y-4">
						<PropertiesPanel />
					</div>
					<div className="card-base p-4 space-y-3">
						<h3 className="text-base font-semibold text-black">GNASS Insights</h3>
						{gnassStatus === 'running' && <p className="text-sm text-midGrey">Training on latest geometry…</p>}
						{gnassStatus === 'error' && (
							<p className="text-sm text-red-500">{gnassError || 'Unable to run GNASS pipeline.'}</p>
						)}
						{gnassStatus === 'ready' && gnassSummary && (
							<ul className="space-y-2">
								{Object.entries(gnassSummary).map(([key, value]) => (
									<li key={key} className="flex items-center justify-between text-sm">
										<span className="text-midGrey capitalize">{key}</span>
										<span className="text-black font-semibold">{value.toFixed(2)}</span>
									</li>
								))}
							</ul>
						)}
						{gnassStatus === 'idle' && <p className="text-sm text-midGrey">Add atoms to train the model.</p>}
					</div>
					<div className="card-base p-4 space-y-3">
						<h3 className="text-base font-semibold text-black">Structure Issues</h3>
						<ul className="space-y-2">
							{issueSummary.warnings.map((warning) => (
								<li key={warning} className="text-sm text-midGrey flex items-start gap-2">
									<span className={`mt-1 h-2 w-2 rounded-full ${issueSummary.hasInstabilityRisk ? 'bg-red-500' : 'bg-emerald-400'}`} />
									<span>{warning}</span>
								</li>
							))}
						</ul>
					</div>
					<div className="card-base p-4 space-y-3">
						<h3 className="text-base font-semibold text-black">Knowledge & Analysis</h3>
						<KABPanel />
					</div>
				</section>
			</div>
		</motion.div>
	);
}


