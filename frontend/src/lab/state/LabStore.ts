/**
 * LabStore - Main Zustand store for Lab state
 */

import { create } from 'zustand';
import { MoleculeStateEngine, type MoleculeState } from '../engines/MoleculeStateEngine';
import { ActionManager } from '../engines/actions/ActionManager';
import { actionRegistry } from '../engines/actions/ActionRegistry';
import { ValidationEngine } from '../engines/validation/ValidationEngine';
import { ToolManager } from '../engines/tools/ToolManager';
import { AtomTool, BondTool, DragTool, EraserTool, SelectTool } from '../engines/tools';
import type { ValidationResult } from '../engines/validation/Validation.types';

interface LabState {
  // Engines
  moleculeEngine: MoleculeStateEngine;
  actionManager: ActionManager;
  validationEngine: ValidationEngine;
  toolManager: ToolManager;
  
  // State
  currentMolecule: MoleculeState | null;
  validationResult: ValidationResult | null;
  selectedAtomId: string | null;
  selectedBondId: string | null;
  currentElement: string;
  
  // Actions
  dispatch: (type: string, payload: any) => void;
  undo: () => void;
  redo: () => void;
  setTool: (toolId: string) => void;
  selectAtom: (atomId: string | null) => void;
  selectBond: (bondId: string | null) => void;
  setCurrentElement: (element: string) => void;
  validate: () => void;
  clear: () => void;
}

export const useLabStore = create<LabState>((set, get) => {
  // Initialize engines
  const moleculeEngine = new MoleculeStateEngine();
  const actionManager = new ActionManager();
  const validationEngine = new ValidationEngine();
  const toolManager = new ToolManager();
  
  // Create selection manager that will be populated after store creation
  let selectionManager: any = null;

  const store = {
    moleculeEngine,
    actionManager,
    validationEngine,
    toolManager,
    currentMolecule: moleculeEngine.getState(),
    validationResult: null,
    selectedAtomId: null,
    selectedBondId: null,
    currentElement: 'C',
    
    dispatch: (type: string, payload: any) => {
      const action = actionRegistry.create(type, payload);
      const currentState = get().currentMolecule || moleculeEngine.getState();
      const newState = actionManager.apply(action, currentState);
      
      // Update molecule engine
      moleculeEngine.loadState(newState);
      
      // Validate
      const validation = validationEngine.validate(newState);
      
      set({
        currentMolecule: newState,
        validationResult: validation,
      });
    },
    
    undo: () => {
      const currentState = get().currentMolecule || moleculeEngine.getState();
      const newState = actionManager.undo(currentState);
      moleculeEngine.loadState(newState);
      const validation = validationEngine.validate(newState);
      
      set({
        currentMolecule: newState,
        validationResult: validation,
      });
    },
    
    redo: () => {
      const currentState = get().currentMolecule || moleculeEngine.getState();
      const newState = actionManager.redo(currentState);
      moleculeEngine.loadState(newState);
      const validation = validationEngine.validate(newState);
      
      set({
        currentMolecule: newState,
        validationResult: validation,
      });
    },
    
    setTool: (toolId: string) => {
      toolManager.setActive(toolId);
      set({});
    },
    
    selectAtom: (atomId: string | null) => {
      set({ selectedAtomId: atomId, selectedBondId: null });
    },
    
    selectBond: (bondId: string | null) => {
      set({ selectedBondId: bondId, selectedAtomId: null });
    },
    
    setCurrentElement: (element: string) => {
      set({ currentElement: element });
    },
    
    validate: () => {
      const state = get().currentMolecule || moleculeEngine.getState();
      const validation = validationEngine.validate(state);
      set({ validationResult: validation });
    },
    
    clear: () => {
      moleculeEngine.clear();
      actionManager.clear();
      set({
        currentMolecule: moleculeEngine.getState(),
        validationResult: null,
        selectedAtomId: null,
        selectedBondId: null,
      });
    },
  };

  // Initialize selection manager and tools after store is created
  selectionManager = {
    selectAtom: (id: string | null) => store.selectAtom(id),
    selectBond: (id: string | null) => store.selectBond(id),
    clear: () => {
      store.selectAtom(null);
      store.selectBond(null);
    },
  };
  
  toolManager.register(SelectTool(moleculeEngine, selectionManager));
  toolManager.register(AtomTool(moleculeEngine));
  toolManager.register(BondTool(moleculeEngine));
  toolManager.register(DragTool(moleculeEngine));
  toolManager.register(EraserTool(moleculeEngine));
  
  // Set default tool
  toolManager.setActive('select');

  return store;
});

