# Lab Architecture - Complete Redesign

This directory contains the complete redesign of the Molecule Lab following a clean, modular architecture.

## Architecture Overview

```
UI (React) → LabStore (Zustand) → MoleculeStateEngine → Actions → Validation
                ↓
            ToolManager → Tools → Actions
```

## Core Components

### 1. Engines (`engines/`)

#### MoleculeStateEngine
- Single source of truth for molecule data
- Manages atoms and bonds as Maps
- Provides methods for adding/removing/updating atoms and bonds
- Serialization to SMILES and SDF (placeholders for backend integration)

#### Action System (`engines/actions/`)
- **Action.types.ts**: Base action interface
- **ActionManager.ts**: Undo/redo stack management
- **ActionRegistry.ts**: Action type registration
- **Actions**: AddAtom, AddBond, DeleteAtom, DeleteBond, MoveAtom

All molecule changes go through actions, ensuring:
- Full undo/redo support
- Action logging
- Serializable history

#### Tool System (`engines/tools/`)
- **ToolManager**: Manages active tool
- **Tools**: SelectTool, AtomTool, BondTool, DragTool, EraserTool
- Tools dispatch actions through the state engine
- Clean separation: tools don't mutate state directly

#### Validation Engine (`engines/validation/`)
- **ValidationEngine**: Main validation orchestrator
- **BasicValenceValidator**: Checks atom valences
- **Sanitizer**: Cleans up molecule state (placeholder)
- Validates after every action

### 2. State Management (`state/`)

#### LabStore (Zustand)
- Central store integrating all engines
- Provides hooks for React components
- Manages:
  - Current molecule state
  - Active tool
  - Selection (atoms/bonds)
  - Validation results
  - Current element for atom tool

### 3. Hooks (`hooks/`)

- **useLab**: Main hook for accessing Lab store
- **useToolMode**: Access tool manager and active tool
- **useSelection**: Access selection state

### 4. UI Components (`components/`)

#### Layout
- **LabPage**: Main page component
- **Toolbar**: Top toolbar with tool buttons and undo/redo
- **LeftSidebar**: Left sidebar with tools and inspector
- **BottomDock**: Bottom panel for validation and logs

#### Editor
- **Editor2D**: 2D canvas-based editor (placeholder - ready for Konva integration)
- **Viewer3D**: 3D viewer placeholder (ready for 3Dmol.js or react-three-fiber)

#### Panels
- **ToolsPanel**: Tool selection panel
- **InspectorPanel**: Atom/bond inspection panel

## Usage

### Basic Example

```tsx
import LabPage from './lab/components/LabPage';

function App() {
  return <LabPage />;
}
```

### Programmatic Usage

```tsx
import { useLab } from './lab/hooks/useLab';

function MyComponent() {
  const { dispatch, undo, redo, moleculeEngine } = useLab();
  
  // Add an atom
  dispatch('addAtom', {
    atomId: crypto.randomUUID(),
    element: 'C',
    position: { x: 100, y: 100 },
  });
  
  // Undo
  undo();
  
  // Get current state
  const atoms = moleculeEngine.getAtoms();
}
```

## File Structure

```
lab/
├── engines/
│   ├── MoleculeStateEngine.ts      # Core molecule state
│   ├── actions/                    # Action system
│   │   ├── Action.types.ts
│   │   ├── ActionManager.ts
│   │   ├── ActionRegistry.ts
│   │   ├── AddAtomAction.ts
│   │   ├── AddBondAction.ts
│   │   ├── DeleteAtomAction.ts
│   │   ├── DeleteBondAction.ts
│   │   ├── MoveAtomAction.ts
│   │   └── index.ts
│   ├── tools/                      # Tool system
│   │   ├── Tool.types.ts
│   │   ├── ToolManager.ts
│   │   ├── AtomTool.ts
│   │   ├── BondTool.ts
│   │   ├── DragTool.ts
│   │   ├── EraserTool.ts
│   │   ├── SelectTool.ts
│   │   └── index.ts
│   └── validation/                 # Validation engine
│       ├── Validation.types.ts
│       ├── ValidationEngine.ts
│       ├── BasicValenceValidator.ts
│       ├── Sanitizer.ts
│       └── index.ts
├── state/
│   └── LabStore.ts                 # Zustand store
├── hooks/                          # React hooks
│   ├── useLab.ts
│   ├── useToolMode.ts
│   └── useSelection.ts
├── components/                     # UI components
│   ├── LabPage.tsx
│   ├── Toolbar.tsx
│   ├── LeftSidebar.tsx
│   ├── Editor2D.tsx
│   ├── Viewer3D.tsx
│   ├── BottomDock.tsx
│   └── panels/
│       ├── ToolsPanel.tsx
│       └── InspectorPanel.tsx
├── index.ts                         # Main exports
└── README.md                        # This file
```

## Next Steps

### Phase 1 Remaining Tasks
1. ✅ Core engines - DONE
2. ✅ Action system - DONE
3. ✅ Tool system - DONE
4. ✅ Validation - DONE
5. ⬜ Full Konva 2D rendering (Editor2D needs proper Konva integration)
6. ⬜ 3D viewer integration (3Dmol.js or react-three-fiber)

### Phase 2: Chemistry Engines
- Molecule serialization (SMILES, MOL, JSON)
- Layout engine (2D coordinate generation)
- Bond angle + hybridization
- 3D geometry generator
- RDKit integration

### Phase 3: Advanced Features
- Spectroscopy engines
- ML model integration
- Export functionality

## Integration with Existing Code

This new architecture is designed to coexist with the existing Lab implementation. To integrate:

1. Import `LabPage` from `./lab/components/LabPage`
2. Add route in your router
3. The new Lab uses its own state management, separate from the existing `moleculeStore`

## Notes

- All actions are reversible (undo/redo)
- Validation runs automatically after each action
- Tools are isolated and don't mutate state directly
- State is centralized in Zustand store
- Ready for Konva.js integration in Editor2D
- Ready for 3Dmol.js or react-three-fiber in Viewer3D

