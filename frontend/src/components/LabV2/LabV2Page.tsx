import { useEffect } from "react";
import { useLocation, useSearchParams } from "react-router-dom";
import TopActionBar from "./layout/TopActionBar";
import LeftToolDock from "./layout/LeftToolDock";
import RightInspector from "./layout/RightInspector";
import BottomTemplateBar from "./layout/BottomTemplateBar";
import LabCanvas from "./LabCanvas";
import Navbar from "../Navbar";
import { useLabStore } from "../../store/labStore";
import type { PublicMolecule, UserMolecule } from "../../lib/userMoleculeStore";
import { getUserMolecule } from "../../lib/userMoleculeStore";
import { getPublicMolecule } from "../../lib/publicMoleculeStore";
import { supabase } from "../../supabase";

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
                return;
            }

            try {
                // If we were passed a molfile directly, prefer that for instant load
                if (state?.molfile) {
                    loadMolecule({
                        id: String(moleculeId),
                        atoms: [],
                        bonds: [],
                        metadata: {
                            name: state.name,
                            smiles: state.smiles,
                            formula: state.formula,
                            molfile: state.molfile,
                            source,
                        },
                    } as any);
                    return;
                }

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
                        loadMolecule({
                            id: userMol.id || String(moleculeId),
                            atoms: [],
                            bonds: [],
                            metadata: {
                                name: userMol.name,
                                smiles: userMol.smiles,
                                formula: userMol.formula,
                                molfile: userMol.molfile,
                                source: "user",
                            },
                        } as any);
                    }
                } else if (source === "public") {
                    const publicMol = await getPublicMolecule(String(moleculeId));
                    if (publicMol) {
                        loadMolecule({
                            id: publicMol.id || String(moleculeId),
                            atoms: [],
                            bonds: [],
                            metadata: {
                                name: publicMol.name,
                                smiles: publicMol.smiles,
                                formula: publicMol.formula,
                                molfile: publicMol.molfile,
                                source: "public",
                            },
                        } as any);
                    }
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
