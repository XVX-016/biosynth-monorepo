import { useState, useRef } from "react";
import {
    MousePointer,
    Circle,
    Minus,
    Move,
    Trash2,
    ZoomIn,
    ZoomOut,
    RotateCcw,
    ChevronRight
} from "lucide-react";
import { useMoleculeStore } from "../../store/moleculeStore";
import type { Element } from "@biosynth/engine";

const ELEMENT_DATA: Record<Element, { color: string; name: string }> = {
    H: { color: "#ffffff", name: "Hydrogen" },
    C: { color: "#909090", name: "Carbon" },
    N: { color: "#3050f8", name: "Nitrogen" },
    O: { color: "#ff0d0d", name: "Oxygen" },
    F: { color: "#90e050", name: "Fluorine" },
    P: { color: "#ff8000", name: "Phosphorus" },
    S: { color: "#ffff30", name: "Sulfur" },
    Cl: { color: "#1ff01f", name: "Chlorine" },
    Br: { color: "#a62929", name: "Bromine" },
};

export function LabToolbar() {
    const tool = useMoleculeStore((state) => state.tool);
    const setTool = useMoleculeStore((state) => state.setTool);
    const setAtomToAdd = useMoleculeStore((state) => state.setAtomToAdd);
    const autoBond = useMoleculeStore((state) => state.autoBond);
    const setAutoBond = useMoleculeStore((state) => state.setAutoBond);
    const currentMolecule = useMoleculeStore((state) => state.currentMolecule);
    const reset = useMoleculeStore((state) => state.reset);

    const [showElementPicker, setShowElementPicker] = useState(false);
    const [selectedElement, setSelectedElement] = useState<Element>("C");

    const commonElements: Element[] = ["C", "H", "O", "N", "S", "P", "F", "Cl", "Br"];

    const handleToolClick = (newTool: "select" | "add-atom" | "bond" | "delete") => {
        if (newTool === "add-atom") {
            setTool("add-atom");
            setAtomToAdd(selectedElement);
            setShowElementPicker(true);
        } else {
            setTool(newTool);
            setShowElementPicker(false);
        }
    };

    const handleElementSelect = (element: Element) => {
        setSelectedElement(element);
        setAtomToAdd(element);
    };

    const tools = [
        {
            id: "select" as const,
            label: "Select",
            icon: MousePointer,
            onClick: () => handleToolClick("select")
        },
        {
            id: "add-atom" as const,
            label: "Atom",
            icon: Circle,
            onClick: () => handleToolClick("add-atom")
        },
        {
            id: "bond" as const,
            label: "Bond",
            icon: Minus,
            onClick: () => handleToolClick("bond")
        },
        {
            id: "delete" as const,
            label: "Delete",
            icon: Trash2,
            onClick: () => handleToolClick("delete")
        },
    ];

    // Zoom functions - will be connected to OrbitControls via ref
    const handleZoomIn = () => {
        const event = new CustomEvent('lab-zoom', { detail: { direction: 'in' } });
        window.dispatchEvent(event);
    };

    const handleZoomOut = () => {
        const event = new CustomEvent('lab-zoom', { detail: { direction: 'out' } });
        window.dispatchEvent(event);
    };

    const handleResetView = () => {
        const event = new CustomEvent('lab-reset-camera');
        window.dispatchEvent(event);
    };

    const handleClearMolecule = () => {
        if (confirm('Clear the entire molecule?')) {
            reset();
        }
    };

    const atomCount = currentMolecule?.atoms.size || 0;
    const bondCount = currentMolecule?.bonds.size || 0;

    return (
        <div className="flex flex-col items-center py-4 gap-2 h-full relative">
            {/* Molecule Stats */}
            {atomCount > 0 && (
                <div className="w-14 px-1 py-2 bg-gray-50 rounded-lg border border-gray-200 mb-2">
                    <div className="text-center">
                        <div className="text-xs font-bold text-gray-700">{atomCount}</div>
                        <div className="text-[9px] text-gray-500">atoms</div>
                    </div>
                    <div className="text-center mt-1">
                        <div className="text-xs font-bold text-gray-700">{bondCount}</div>
                        <div className="text-[9px] text-gray-500">bonds</div>
                    </div>
                </div>
            )}

            {/* Tool Buttons */}
            <div className="flex flex-col gap-2">
                {tools.map((toolItem) => {
                    const Icon = toolItem.icon;
                    const isActive = tool === toolItem.id;

                    return (
                        <button
                            key={toolItem.id}
                            onClick={toolItem.onClick}
                            className={`
                relative flex flex-col items-center justify-center w-14 h-14 rounded-lg
                transition-all duration-200
                ${isActive
                                    ? "bg-black text-white shadow-lg"
                                    : "bg-gray-100 text-gray-700 hover:bg-gray-200"
                                }
              `}
                            title={toolItem.label}
                        >
                            <Icon size={20} />
                            <span className="text-[10px] mt-1">{toolItem.label}</span>
                            {toolItem.id === "add-atom" && isActive && (
                                <div className="absolute -right-1 -top-1 w-3 h-3 bg-blue-500 rounded-full animate-pulse" />
                            )}
                        </button>
                    );
                })}
            </div>

            {/* Enhanced Element Picker */}
            {showElementPicker && tool === "add-atom" && (
                <div className="absolute left-20 top-4 bg-white rounded-xl shadow-2xl border border-gray-200 z-50 overflow-hidden">
                    {/* Header */}
                    <div className="bg-gradient-to-r from-gray-900 to-gray-800 px-4 py-3 text-white">
                        <div className="text-sm font-bold">Element Selector</div>
                        <div className="text-xs opacity-75 mt-0.5">Click to add atoms</div>
                    </div>

                    {/* Current Selection */}
                    <div className="px-4 py-3 bg-gray-50 border-b border-gray-200">
                        <div className="text-xs text-gray-600 mb-2">Selected Element:</div>
                        <div className="flex items-center gap-3">
                            <div
                                className="w-10 h-10 rounded-lg flex items-center justify-center font-bold text-white shadow-md"
                                style={{ backgroundColor: ELEMENT_DATA[selectedElement].color === "#ffffff" ? "#e0e0e0" : ELEMENT_DATA[selectedElement].color }}
                            >
                                {selectedElement}
                            </div>
                            <div>
                                <div className="font-semibold text-gray-800">{ELEMENT_DATA[selectedElement].name}</div>
                                <div className="text-xs text-gray-500">Click canvas to place</div>
                            </div>
                        </div>
                    </div>

                    {/* Element Grid */}
                    <div className="p-4">
                        <div className="grid grid-cols-3 gap-2">
                            {commonElements.map((element) => {
                                const data = ELEMENT_DATA[element];
                                const isSelected = selectedElement === element;

                                return (
                                    <button
                                        key={element}
                                        onClick={() => handleElementSelect(element)}
                                        className={`
                      group relative px-3 py-3 rounded-lg font-mono font-bold text-sm
                      transition-all duration-200
                      ${isSelected
                                                ? "ring-2 ring-blue-500 ring-offset-2 scale-105"
                                                : "hover:scale-105 hover:shadow-md"
                                            }
                    `}
                                        style={{
                                            backgroundColor: data.color === "#ffffff" ? "#f0f0f0" : data.color,
                                            color: ["H", "F", "Cl", "S"].includes(element) ? "#000" : "#fff"
                                        }}
                                        title={data.name}
                                    >
                                        <div className="text-center">{element}</div>
                                        {isSelected && (
                                            <div className="absolute -top-1 -right-1 w-2 h-2 bg-blue-500 rounded-full" />
                                        )}
                                    </button>
                                );
                            })}
                        </div>
                    </div>

                    {/* Auto-bond Settings */}
                    <div className="px-4 py-3 bg-gray-50 border-t border-gray-200">
                        <label className="flex items-center gap-2 cursor-pointer group">
                            <div className="relative">
                                <input
                                    type="checkbox"
                                    checked={autoBond}
                                    onChange={(e) => setAutoBond(e.target.checked)}
                                    className="w-4 h-4 rounded border-gray-300 text-blue-600 focus:ring-blue-500"
                                />
                            </div>
                            <div className="flex-1">
                                <div className="text-sm font-medium text-gray-800">Auto-Bond</div>
                                <div className="text-xs text-gray-500">ML-powered bond formation</div>
                            </div>
                            <div className={`w-2 h-2 rounded-full ${autoBond ? 'bg-green-500' : 'bg-gray-300'}`} />
                        </label>
                    </div>
                </div>
            )}

            {/* Spacer */}
            <div className="flex-1" />

            {/* Action Buttons */}
            {atomCount > 0 && (
                <button
                    onClick={handleClearMolecule}
                    className="w-14 p-2 rounded-lg bg-red-50 hover:bg-red-100 text-red-600 transition-colors border border-red-200"
                    title="Clear All"
                >
                    <Trash2 size={16} className="mx-auto" />
                    <span className="text-[9px] mt-1 block">Clear</span>
                </button>
            )}

            {/* View Controls */}
            <div className="flex flex-col gap-2 border-t border-gray-200 pt-2">
                <button
                    onClick={handleZoomIn}
                    className="p-3 rounded-lg bg-gray-100 hover:bg-gray-200 transition-colors"
                    title="Zoom In"
                >
                    <ZoomIn size={20} className="text-gray-700" />
                </button>
                <button
                    onClick={handleZoomOut}
                    className="p-3 rounded-lg bg-gray-100 hover:bg-gray-200 transition-colors"
                    title="Zoom Out"
                >
                    <ZoomOut size={20} className="text-gray-700" />
                </button>
                <button
                    onClick={handleResetView}
                    className="p-3 rounded-lg bg-gray-100 hover:bg-gray-200 transition-colors"
                    title="Reset View"
                >
                    <RotateCcw size={20} className="text-gray-700" />
                </button>
            </div>
        </div>
    );
}
