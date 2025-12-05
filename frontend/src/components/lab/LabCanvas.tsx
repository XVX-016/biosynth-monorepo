import { useRef, useEffect } from "react";
import { Canvas, useThree } from "@react-three/fiber";
import { OrbitControls, Grid, Line } from "@react-three/drei";
import { useMoleculeStore } from "../../store/moleculeStore";
import { addAtom, moleculeToRenderable } from "../../lib/engineAdapter";
import { Vector3, Vector2, Raycaster } from "three";
import AtomMesh from "../r3f/AtomMesh";
import BondMesh from "../r3f/BondMesh";
import type { OrbitControls as OrbitControlsType } from "three-stdlib";

// Camera controller for toolbar buttons
function CameraController() {
    const { camera, controls } = useThree();
    const controlsRef = useRef<OrbitControlsType | null>(null);

    useEffect(() => {
        if (controls) {
            controlsRef.current = controls as OrbitControlsType;
        }
    }, [controls]);

    useEffect(() => {
        const handleZoom = (event: Event) => {
            const customEvent = event as CustomEvent<{ direction: 'in' | 'out' }>;
            const { direction } = customEvent.detail;

            if (controlsRef.current) {
                const zoomFactor = direction === 'in' ? 0.8 : 1.2;
                camera.position.multiplyScalar(zoomFactor);
                controlsRef.current.update();
            }
        };

        const handleReset = () => {
            if (controlsRef.current) {
                camera.position.set(15, 12, 15);
                controlsRef.current.target.set(0, 0, 0);
                controlsRef.current.update();
            }
        };

        window.addEventListener('lab-zoom', handleZoom);
        window.addEventListener('lab-reset-camera', handleReset);

        return () => {
            window.removeEventListener('lab-zoom', handleZoom);
            window.removeEventListener('lab-reset-camera', handleReset);
        };
    }, [camera]);

    return null;
}

// Reference Axes and Planes
function ReferenceAxes() {
    const gridSize = 50;
    const axisColor = "#a0a0a0";

    return (
        <group>
            {/* Vertical reference lines at edges */}
            {/* Left edge */}
            <Line
                points={[[-gridSize, -5, 0], [-gridSize, 5, 0]]}
                color={axisColor}
                lineWidth={1}
                dashed={false}
            />
            {/* Right edge */}
            <Line
                points={[[gridSize, -5, 0], [gridSize, 5, 0]]}
                color={axisColor}
                lineWidth={1}
                dashed={false}
            />
            {/* Front edge */}
            <Line
                points={[[0, -5, -gridSize], [0, 5, -gridSize]]}
                color={axisColor}
                lineWidth={1}
                dashed={false}
            />
            {/* Back edge */}
            <Line
                points={[[0, -5, gridSize], [0, 5, gridSize]]}
                color={axisColor}
                lineWidth={1}
                dashed={false}
            />
        </group>
    );
}

// Double-sided infinite grid
function DoubleSidedGrid() {
    return (
        <>
            {/* Top-facing grid */}
            <Grid
                args={[300, 300]}
                cellSize={1}
                cellThickness={0.5}
                cellColor="#e8e8e8"
                sectionSize={5}
                sectionThickness={0.8}
                sectionColor="#d0d0d0"
                fadeDistance={150}
                fadeStrength={1.5}
                followCamera={false}
                infiniteGrid={true}
                position={[0, 0, 0]}
                rotation={[0, 0, 0]}
            />

            {/* Bottom-facing grid (visible from underneath) */}
            <Grid
                args={[300, 300]}
                cellSize={1}
                cellThickness={0.5}
                cellColor="#e8e8e8"
                sectionSize={5}
                sectionThickness={0.8}
                sectionColor="#d0d0d0"
                fadeDistance={150}
                fadeStrength={1.5}
                followCamera={false}
                infiniteGrid={true}
                position={[0, 0, 0]}
                rotation={[Math.PI, 0, 0]}
            />
        </>
    );
}

// Click handler for adding atoms with auto-bond validation
function CanvasClickHandler() {
    const { camera, gl } = useThree();
    const tool = useMoleculeStore((state) => state.tool);
    const atomToAdd = useMoleculeStore((state) => state.atomToAdd);
    const autoBond = useMoleculeStore((state) => state.autoBond);
    const runValidation = useMoleculeStore((state) => state.runValidation);

    useEffect(() => {
        const canvas = gl.domElement;

        const handleClick = async (event: MouseEvent) => {
            if (tool !== "add-atom" || !atomToAdd) return;

            // Get click position in normalized device coordinates
            const rect = canvas.getBoundingClientRect();
            const x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
            const y = -((event.clientY - rect.top) / rect.height) * 2 + 1;

            // Raycast to the ground plane (y = 0)
            const raycaster = new Raycaster();
            const mousePos = new Vector2(x, y);
            raycaster.setFromCamera(mousePos, camera);

            // Intersect with ground plane at y=0
            const planeY = 0;
            const direction = raycaster.ray.direction;
            const origin = raycaster.ray.origin;

            if (Math.abs(direction.y) < 0.0001) return; // Ray parallel to plane

            const t = (planeY - origin.y) / direction.y;
            if (t < 0) return; // Behind camera

            const intersectPoint = new Vector3(
                origin.x + direction.x * t,
                planeY,
                origin.z + direction.z * t
            );

            console.log('[LabCanvas] Adding atom:', atomToAdd, 'at position:', intersectPoint, 'autoBond:', autoBond);

            // Add atom at intersection point
            addAtom(atomToAdd, [intersectPoint.x, intersectPoint.y, intersectPoint.z]);

            // Run validation after adding atom
            if (autoBond) {
                setTimeout(() => {
                    runValidation();
                }, 200);
            }
        };

        canvas.addEventListener("click", handleClick);
        return () => canvas.removeEventListener("click", handleClick);
    }, [tool, atomToAdd, autoBond, camera, gl, runValidation]);

    return null;
}

// Scene content with molecules and validation highlighting
function SceneContent() {
    const currentMolecule = useMoleculeStore((state) => state.currentMolecule);
    const renderMode = useMoleculeStore((state) => state.renderMode);
    const colorScheme = useMoleculeStore((state) => state.colorScheme);
    const highlightedAtoms = useMoleculeStore((state) => state.highlightedAtoms);
    const validation = useMoleculeStore((state) => state.validation);

    const { atoms, bonds } = moleculeToRenderable(currentMolecule);

    // Get atoms with validation issues
    const atomsWithIssues = new Set<string>();
    if (validation.result?.issues) {
        validation.result.issues.forEach(issue => {
            issue.atoms?.forEach(atomId => atomsWithIssues.add(atomId));
        });
    }

    return (
        <>
            {atoms.map((atom) => {
                const hasIssue = atomsWithIssues.has(atom.id);
                const isHighlighted = highlightedAtoms.includes(atom.id) || hasIssue;

                return (
                    <AtomMesh
                        key={atom.id}
                        id={atom.id}
                        position={atom.position}
                        element={atom.element as any}
                        renderMode={renderMode}
                        colorScheme={colorScheme}
                        highlighted={isHighlighted}
                    />
                );
            })}
            {bonds.map((bond) => (
                <BondMesh
                    key={bond.id}
                    id={bond.id}
                    from={bond.from}
                    to={bond.to}
                    order={bond.order}
                    renderMode={renderMode}
                />
            ))}
        </>
    );
}

// Validation Status Indicator
function ValidationIndicator() {
    const validation = useMoleculeStore((state) => state.validation);
    const currentMolecule = useMoleculeStore((state) => state.currentMolecule);
    const autoBond = useMoleculeStore((state) => state.autoBond);

    if (!currentMolecule || currentMolecule.atoms.size === 0) return null;

    const hasIssues = validation.result && validation.result.issues.length > 0;
    const score = validation.result?.score ?? 100;

    return (
        <div className="absolute top-4 right-4 bg-white/90 backdrop-blur-sm rounded-lg shadow-lg border border-gray-200 p-3 min-w-[200px]">
            <div className="flex items-center justify-between mb-2">
                <span className="text-xs font-semibold text-gray-700">Structure Validation</span>
                <div className="flex items-center gap-2">
                    {autoBond && (
                        <span className="text-[10px] px-2 py-0.5 bg-blue-100 text-blue-700 rounded">Auto-Bond ON</span>
                    )}
                    <div className={`w-2 h-2 rounded-full ${score >= 80 ? 'bg-green-500' : score >= 50 ? 'bg-yellow-500' : 'bg-red-500'}`} />
                </div>
            </div>

            <div className="mb-2">
                <div className="flex justify-between text-xs mb-1">
                    <span className="text-gray-600">Quality Score</span>
                    <span className="font-bold text-gray-800">{score}/100</span>
                </div>
                <div className="w-full bg-gray-200 rounded-full h-1.5">
                    <div
                        className={`h-1.5 rounded-full transition-all ${score >= 80 ? 'bg-green-500' : score >= 50 ? 'bg-yellow-500' : 'bg-red-500'}`}
                        style={{ width: `${score}%` }}
                    />
                </div>
            </div>

            {hasIssues && (
                <div className="text-xs">
                    <div className="text-red-600 font-medium">
                        {validation.result!.issues.length} issue{validation.result!.issues.length !== 1 ? 's' : ''} detected
                    </div>
                    <div className="text-gray-500 mt-1">
                        {validation.result!.suggestions[0]}
                    </div>
                </div>
            )}

            {!hasIssues && score === 100 && (
                <div className="text-xs text-green-600 font-medium">
                    ✓ Structure is valid
                </div>
            )}
        </div>
    );
}

export function LabCanvas() {
    return (
        <div className="w-full h-full relative">
            <Canvas
                camera={{ position: [15, 12, 15], fov: 50 }}
                style={{ width: "100%", height: "100%", background: "#ffffff" }}
            >
                {/* Minimal Lighting - Very subtle */}
                <ambientLight intensity={0.95} />
                <directionalLight intensity={0.15} position={[10, 20, 10]} />

                {/* Double-sided Grid - Visible from both sides */}
                <DoubleSidedGrid />

                {/* Reference axes at edges */}
                <ReferenceAxes />

                {/* Molecule rendering with validation */}
                <SceneContent />

                {/* Click handler for atom placement */}
                <CanvasClickHandler />

                {/* Camera controller */}
                <CameraController />

                {/* Interactive Controls - TRUE 3D (can go underneath) */}
                <OrbitControls
                    enableDamping
                    dampingFactor={0.05}
                    enableZoom={true}
                    enablePan={true}
                    enableRotate={true}
                    minDistance={3}
                    maxDistance={100}
                    target={[0, 0, 0]}
                />
            </Canvas>

            {/* Validation Indicator Overlay */}
            <ValidationIndicator />
        </div>
    );
}
