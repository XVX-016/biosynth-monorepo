import React, { useState, useEffect } from 'react';
import { useSearchParams } from 'react-router-dom';
import MoleculeViewer from '../components/MoleculeViewer';
import ToolPanel from '../components/ToolPanel';
import PropertiesPanel from '../components/PropertiesPanel';
import { TemplatePanel } from '../components/TemplatePanel';
import Button from '../components/ui/Button';
import { useMoleculeStore } from '../store/moleculeStore';
import { saveUserMolecule, getUserMolecule } from '../lib/userMoleculeStore';
import { getPublicMolecule } from '../lib/publicMoleculeStore';
import { moleculeToJSON, getCanvasThumbnail, moleculeFromJSON } from '../lib/engineAdapter';
import { MoleculeSerializer, MoleculeGraph } from '@biosynth/engine';
import { supabase } from '../supabase';
import { convertSMILESToMolfile, generateThumbnailBase64 } from '../lib/api';

export default function Lab() {
	const [searchParams] = useSearchParams();
	const molecule = useMoleculeStore((s) => s.currentMolecule);
	const setMolecule = useMoleculeStore((s) => s.setMolecule);
	const backendPredictions = useMoleculeStore((s) => s.backendPredictions);
	const [saving, setSaving] = useState(false);
	const [userId, setUserId] = useState<string | null>(null);
	const [loadingMolecule, setLoadingMolecule] = useState(false);

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
		<div className="lab-layout">
			{loadingMolecule && (
				<div className="fixed inset-0 bg-black/50 flex items-center justify-center z-50">
					<div className="bg-white rounded-lg p-6">
						<div className="text-center">Loading molecule...</div>
					</div>
				</div>
			)}
			<TemplatePanel />
			<div className="grid grid-cols-12 gap-4 flex-1">
				<div className="col-span-12 lg:col-span-3">
					<div className="bg-white rounded-xl shadow-neon border border-lightGrey p-4">
						<ToolPanel />
					</div>
				</div>
				<div className="col-span-12 lg:col-span-6">
					<div className="relative rounded-xl shadow-neon border border-lightGrey p-2 h-[70vh] lg:h-[78vh] bg-offwhite">
						<MoleculeViewer />
						<div className="absolute bottom-3 right-3">
							<Button onClick={handleSave} disabled={!molecule || saving || !userId}>
								{saving ? 'Saving...' : userId ? 'Save to Library' : 'Sign in to Save'}
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


