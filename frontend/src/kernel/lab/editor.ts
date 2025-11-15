import { MoleculeGraph } from '@biosynth/engine';
import { UndoStack } from '@biosynth/engine';

export type ManipulationMode = 'drag' | 'rotate' | 'bond' | 'select' | 'measure';

export interface SelectionState {
  selectedAtomIds: Set<string>;
  selectedBondIds: Set<string>;
  hoveredAtomId: string | null;
  hoveredBondId: string | null;
}

export interface Measurement {
  id: string;
  type: 'distance' | 'angle' | 'dihedral';
  atomIds: string[];
  value: number;
  unit: string;
}

export interface ClipboardItem {
  type: 'atom' | 'fragment';
  data: MoleculeGraph | { element: string; position: [number, number, number] };
}

class LabEditor {
  private selectionState: SelectionState = {
    selectedAtomIds: new Set(),
    selectedBondIds: new Set(),
    hoveredAtomId: null,
    hoveredBondId: null,
  };
  
  private manipulationMode: ManipulationMode = 'select';
  private undoStack: UndoStack<MoleculeGraph>;
  private measurements: Measurement[] = [];
  private clipboard: ClipboardItem | null = null;
  
  constructor() {
    this.undoStack = new UndoStack<MoleculeGraph>();
  }
  
  // Selection management
  selectAtom(atomId: string, multiSelect: boolean = false): void {
    if (multiSelect) {
      if (this.selectionState.selectedAtomIds.has(atomId)) {
        this.selectionState.selectedAtomIds.delete(atomId);
      } else {
        this.selectionState.selectedAtomIds.add(atomId);
      }
    } else {
      this.selectionState.selectedAtomIds.clear();
      this.selectionState.selectedAtomIds.add(atomId);
    }
  }
  
  selectBond(bondId: string, multiSelect: boolean = false): void {
    if (multiSelect) {
      if (this.selectionState.selectedBondIds.has(bondId)) {
        this.selectionState.selectedBondIds.delete(bondId);
      } else {
        this.selectionState.selectedBondIds.add(bondId);
      }
    } else {
      this.selectionState.selectedBondIds.clear();
      this.selectionState.selectedBondIds.add(bondId);
    }
  }
  
  clearSelection(): void {
    this.selectionState.selectedAtomIds.clear();
    this.selectionState.selectedBondIds.clear();
  }
  
  setHoveredAtom(atomId: string | null): void {
    this.selectionState.hoveredAtomId = atomId;
  }
  
  setHoveredBond(bondId: string | null): void {
    this.selectionState.hoveredBondId = bondId;
  }
  
  getSelectionState(): SelectionState {
    return {
      selectedAtomIds: new Set(this.selectionState.selectedAtomIds),
      selectedBondIds: new Set(this.selectionState.selectedBondIds),
      hoveredAtomId: this.selectionState.hoveredAtomId,
      hoveredBondId: this.selectionState.hoveredBondId,
    };
  }
  
  // Manipulation modes
  setManipulationMode(mode: ManipulationMode): void {
    this.manipulationMode = mode;
  }
  
  getManipulationMode(): ManipulationMode {
    return this.manipulationMode;
  }
  
  // Undo/Redo
  pushState(molecule: MoleculeGraph): void {
    this.undoStack.push(molecule.clone());
  }
  
  undo(): MoleculeGraph | null {
    const state = this.undoStack.undo();
    return state ? state.clone() : null;
  }
  
  redo(): MoleculeGraph | null {
    const state = this.undoStack.redo();
    return state ? state.clone() : null;
  }
  
  canUndo(): boolean {
    return this.undoStack.canUndo();
  }
  
  canRedo(): boolean {
    return this.undoStack.canRedo();
  }
  
  clearHistory(): void {
    this.undoStack = new UndoStack<MoleculeGraph>();
  }
  
  // Measurement tools
  addMeasurement(measurement: Measurement): void {
    this.measurements.push(measurement);
  }
  
  removeMeasurement(id: string): void {
    this.measurements = this.measurements.filter((m) => m.id !== id);
  }
  
  clearMeasurements(): void {
    this.measurements = [];
  }
  
  getMeasurements(): Measurement[] {
    return [...this.measurements];
  }
  
  /**
   * Calculate distance between two atoms
   */
  measureDistance(atomId1: string, atomId2: string, molecule: MoleculeGraph): number | null {
    const atom1 = molecule.atoms.get(atomId1);
    const atom2 = molecule.atoms.get(atomId2);
    
    if (!atom1 || !atom2) return null;
    
    const dx = atom1.position[0] - atom2.position[0];
    const dy = atom1.position[1] - atom2.position[1];
    const dz = atom1.position[2] - atom2.position[2];
    
    return Math.sqrt(dx * dx + dy * dy + dz * dz);
  }
  
  /**
   * Calculate angle between three atoms (atom2 is the vertex)
   */
  measureAngle(atomId1: string, atomId2: string, atomId3: string, molecule: MoleculeGraph): number | null {
    const atom1 = molecule.atoms.get(atomId1);
    const atom2 = molecule.atoms.get(atomId2);
    const atom3 = molecule.atoms.get(atomId3);
    
    if (!atom1 || !atom2 || !atom3) return null;
    
    // Vector from atom2 to atom1
    const v1 = [
      atom1.position[0] - atom2.position[0],
      atom1.position[1] - atom2.position[1],
      atom1.position[2] - atom2.position[2],
    ];
    
    // Vector from atom2 to atom3
    const v2 = [
      atom3.position[0] - atom2.position[0],
      atom3.position[1] - atom2.position[1],
      atom3.position[2] - atom2.position[2],
    ];
    
    // Dot product
    const dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    
    // Magnitudes
    const mag1 = Math.sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
    const mag2 = Math.sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
    
    if (mag1 === 0 || mag2 === 0) return null;
    
    // Angle in radians, convert to degrees
    const cosAngle = dot / (mag1 * mag2);
    const angle = Math.acos(Math.max(-1, Math.min(1, cosAngle)));
    
    return (angle * 180) / Math.PI;
  }
  
  // Clipboard
  copyToClipboard(item: ClipboardItem): void {
    this.clipboard = item;
  }
  
  getClipboard(): ClipboardItem | null {
    return this.clipboard;
  }
  
  clearClipboard(): void {
    this.clipboard = null;
  }
  
  /**
   * Copy selected atoms/fragment to clipboard
   */
  copySelection(molecule: MoleculeGraph): void {
    if (this.selectionState.selectedAtomIds.size === 0) return;
    
    // Create fragment from selected atoms
    const fragment = new MoleculeGraph();
    const idMap = new Map<string, string>();
    
    // Add selected atoms
    for (const atomId of this.selectionState.selectedAtomIds) {
      const atom = molecule.atoms.get(atomId);
      if (atom) {
        const newId = fragment.addAtom({
          element: atom.element,
          position: atom.position,
        });
        idMap.set(atomId, newId);
      }
    }
    
    // Add bonds between selected atoms
    for (const bond of molecule.bonds.values()) {
      if (
        this.selectionState.selectedAtomIds.has(bond.a1) &&
        this.selectionState.selectedAtomIds.has(bond.a2)
      ) {
        const newA1 = idMap.get(bond.a1);
        const newA2 = idMap.get(bond.a2);
        if (newA1 && newA2) {
          fragment.addBond(newA1, newA2, bond.order);
        }
      }
    }
    
    this.copyToClipboard({
      type: 'fragment',
      data: fragment,
    });
  }
  
  /**
   * Paste from clipboard
   */
  pasteFromClipboard(
    molecule: MoleculeGraph,
    position: [number, number, number] = [0, 0, 0]
  ): MoleculeGraph | null {
    if (!this.clipboard) return null;
    
    if (this.clipboard.type === 'fragment' && this.clipboard.data instanceof MoleculeGraph) {
      const fragment = this.clipboard.data;
      const merged = molecule.clone();
      const idMap = new Map<string, string>();
      
      // Add atoms with offset
      for (const atom of fragment.atoms.values()) {
        const newPosition: [number, number, number] = [
          atom.position[0] + position[0],
          atom.position[1] + position[1],
          atom.position[2] + position[2],
        ];
        const newId = merged.addAtom({
          element: atom.element,
          position: newPosition,
        });
        idMap.set(atom.id, newId);
      }
      
      // Add bonds
      for (const bond of fragment.bonds.values()) {
        const newA1 = idMap.get(bond.a1);
        const newA2 = idMap.get(bond.a2);
        if (newA1 && newA2) {
          merged.addBond(newA1, newA2, bond.order);
        }
      }
      
      return merged;
    }
    
    return null;
  }
  
  // Reset
  reset(): void {
    this.selectionState = {
      selectedAtomIds: new Set(),
      selectedBondIds: new Set(),
      hoveredAtomId: null,
      hoveredBondId: null,
    };
    this.manipulationMode = 'select';
    this.clearHistory();
    this.clearMeasurements();
    this.clearClipboard();
  }
}

// Singleton instance
export const labEditor = new LabEditor();

