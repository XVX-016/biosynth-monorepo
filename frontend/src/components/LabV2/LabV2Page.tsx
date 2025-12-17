import { useEffect } from "react";
import { useLocation, useSearchParams } from "react-router-dom";
import TopActionBar from "./layout/TopActionBar";
import LeftToolDock from "./layout/LeftToolDock";
import RightInspector from "./layout/RightInspector";
import BottomTemplateBar from "./layout/BottomTemplateBar";
import LabCanvas from "./LabCanvas";
import Navbar from "../Navbar";
import { useLabStore } from "../../store/labStore";
import { getUserMolecule } from "../../lib/userMoleculeStore";
import { getPublicMolecule } from "../../lib/publicMoleculeStore";
import { supabase } from "../../supabase";
import { parseMoleculeFromSupabase } from "../../utils/moleculeParser";

export default function LabV2Page() {
    const location = useLocation();
    const [searchParams] = useSearchParams();
    const loadMolecule = useLabStore(s => s.loadMolecule);

    useEffect(() => {
        const loadFromNavigation = async () => {
            const state = location.state as
                | undefined
                | {
                    source?: "user" | "public";
                    moleculeId?: string | number;
                    molfile?: string | null;
                    name?: string;
                    smiles?: string | null;
                    formula?: string | null;
                };

            let source = state?.source as "user" | "public" | undefined;
            let moleculeId = state?.moleculeId;

            // Fallback to query params if state is not present (e.g. hard reload)
            if (!source) {
                const qpSource = searchParams.get("source");
                if (qpSource === "user" || qpSource === "public") {
                    source = qpSource;
                }
            }
            if (!moleculeId) {
                const qpId = searchParams.get("id");
                if (qpId) {
                    moleculeId = qpId;
                }
            }

            if (!source || !moleculeId) {
                // eslint-disable-next-line no-console
                console.log('[LabV2Page] Missing source or moleculeId:', { source, moleculeId });
                return;
            }

            // eslint-disable-next-line no-console
            console.log('[LabV2Page] Fetching molecule:', { source, moleculeId });

            try {
                let moleculeData: {
                    id?: string | number;
                    name?: string;
                    smiles?: string | null;
                    formula?: string | null;
                    molfile?: string | null;
                    json_graph?: string | null;
                    source?: "user" | "public";
                } | null = null;

                // Fetch molecule from Supabase
                if (source === "user") {
                    if (!supabase) {
                        return;
                    }
                    const { data: sessionData } = await supabase.auth.getSession();
                    const userId = sessionData.session?.user?.id;
                    if (!userId) {
                        return;
                    }
                    const userMol = await getUserMolecule(userId, String(moleculeId));
                    if (userMol) {
                        moleculeData = {
                            id: userMol.id || moleculeId,
                            name: userMol.name,
                            smiles: userMol.smiles,
                            formula: userMol.formula,
                            molfile: userMol.molfile,
                            json_graph: (userMol as any).json_graph,
                            source: "user",
                        };
                        // eslint-disable-next-line no-console
                        console.log('[LabV2Page] Fetched user molecule data:', {
                            id: moleculeData.id,
                            name: moleculeData.name,
                            hasMolfile: !!moleculeData.molfile,
                            hasJsonGraph: !!moleculeData.json_graph,
                        });
                    } else {
                        // eslint-disable-next-line no-console
                        console.warn('[LabV2Page] User molecule not found:', { userId, moleculeId });
                    }
                } else if (source === "public") {
                    const publicMol = await getPublicMolecule(String(moleculeId));
                    if (publicMol) {
                        moleculeData = {
                            id: publicMol.id || moleculeId,
                            name: publicMol.name,
                            smiles: publicMol.smiles,
                            formula: publicMol.formula,
                            molfile: publicMol.molfile,
                            json_graph: (publicMol as any).json_graph,
                            source: "public",
                        };
                        // eslint-disable-next-line no-console
                        console.log('[LabV2Page] Fetched public molecule data:', {
                            id: moleculeData.id,
                            name: moleculeData.name,
                            hasMolfile: !!moleculeData.molfile,
                            hasJsonGraph: !!moleculeData.json_graph,
                        });
                    } else {
                        // eslint-disable-next-line no-console
                        console.warn('[LabV2Page] Public molecule not found:', { moleculeId });
                    }
                }

                // If we have molecule data, parse and load it
                if (moleculeData) {
                    // If molecule has SMILES but no molfile/json_graph, try to generate from SMILES
                    if (!moleculeData.molfile && !moleculeData.json_graph && moleculeData.smiles && moleculeData.smiles.trim()) {
                        // eslint-disable-next-line no-console
                        console.log('[LabV2Page] Generating molfile from SMILES:', moleculeData.smiles);
                        try {
                            // Try backend API first (better quality 3D coordinates)
                            const { convertSMILESToMolfile } = await import('../../lib/api');
                            const result = await convertSMILESToMolfile(moleculeData.smiles.trim());
                            if (result.molfile) {
                                moleculeData.molfile = result.molfile;
                                // eslint-disable-next-line no-console
                                console.log('[LabV2Page] Successfully generated molfile from SMILES via backend API');
                                // Optionally save to Supabase (async, don't block)
                                // This could be done in background to persist for future loads
                            }
                        } catch (e) {
                            // eslint-disable-next-line no-console
                            console.warn('[LabV2Page] Backend SMILES conversion failed, will use frontend parser fallback:', e);
                            // Frontend parser fallback will be handled in parseMoleculeFromSupabase
                        }
                    }

                    // eslint-disable-next-line no-console
                    console.log('[LabV2Page] Parsing molecule from Supabase data');
                    const molecule = parseMoleculeFromSupabase(moleculeData);
                    // eslint-disable-next-line no-console
                    console.log('[LabV2Page] Parsed molecule:', {
                        id: molecule.id,
                        atomCount: molecule.atoms.length,
                        bondCount: molecule.bonds.length,
                        name: molecule.metadata?.name,
                    });
                    // eslint-disable-next-line no-console
                    console.log('[LabV2Page] Loading molecule into store');
                    loadMolecule(molecule);
                } else if (state?.molfile) {
                    // Fallback: if we have molfile in state but no Supabase data, parse it directly
                    // eslint-disable-next-line no-console
                    console.log('[LabV2Page] Using fallback molfile from state');
                    const molecule = parseMoleculeFromSupabase({
                        id: moleculeId,
                        name: state.name,
                        smiles: state.smiles,
                        formula: state.formula,
                        molfile: state.molfile,
                        source,
                    });
                    // eslint-disable-next-line no-console
                    console.log('[LabV2Page] Parsed fallback molecule:', {
                        id: molecule.id,
                        atomCount: molecule.atoms.length,
                        bondCount: molecule.bonds.length,
                    });
                    loadMolecule(molecule);
                } else {
                    // eslint-disable-next-line no-console
                    console.warn('[LabV2Page] No molecule data available to load');
                }
            } catch (e) {
                // eslint-disable-next-line no-console
                console.error("[LabV2Page] Failed to load molecule for Lab", e);
            }
        };

        void loadFromNavigation();
    }, [location.state, searchParams, loadMolecule]);

    return (
        <div className="h-screen w-full grid grid-rows-[64px_64px_1fr_72px] overflow-hidden bg-white">
            {/* Row 1: Global Navbar */}
            <div className="z-30 border-b border-gray-100 bg-white">
                <Navbar />
            </div>

            {/* Row 2: Top Toolbar */}
            <div className="z-20 border-b border-gray-100 bg-white">
                <TopActionBar />
            </div>

            {/* Row 3: Body (Palette + Canvas) */}
            <div className="grid grid-cols-[72px_1fr] min-h-0 relative z-0">
                {/* Col 1: Palette */}
                <div className="h-full z-10">
                    <LeftToolDock />
                </div>

                {/* Col 2: Canvas */}
                <div className="relative w-full h-full bg-[#f8f9fa] overflow-hidden">
                    {/* The Sacred Canvas */}
                    <div className="absolute inset-0 z-0">
                        <LabCanvas />
                    </div>

                    {/* Floating Inspector (Overlay) */}
                    {/* We keep this here to slide over canvas without resizing it, 
                        preserving the 'Maximum uninterrupted canvas' goal until interaction */}
                    <RightInspector />
                </div>
            </div>

            {/* Row 4: Bottom Templates */}
            <div className="z-20 border-t border-gray-100 bg-white">
                <BottomTemplateBar />
            </div>
        </div>
    );
}
