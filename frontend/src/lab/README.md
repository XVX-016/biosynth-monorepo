# Molecule Lab Architecture

This directory contains the complete Molecule Lab system with clean, modular architecture.

## Architecture Overview

```
UI (React) → MoleculeStateEngine → GeometryEngine (RDKit) → AI Models
```

### Core Engines

1. **MoleculeStateEngine** (`engines/MoleculeStateEngine.ts`)
   - Single source of truth for molecule structure
   - Manages atoms and bonds
   - Handles serialization (SMILES, SDF, JSON)
   - No rendering logic - pure state management

2. **ToolController** (`engines/ToolController.ts`)
   - Manages active tool mode (cursor, addAtom, addBond, erase, select)
   - Routes mouse events to appropriate handlers
   - Handles tool-specific state (bond preview, selection)

3. **CommandStack** (`engines/CommandStack.ts`)
   - Implements undo/redo pattern
   - Commands: AddAtom, RemoveAtom, AddBond, MoveAtom
   - Reliable history management

4. **ValidationEngine** (`engines/ValidationEngine.ts`)
   - Local quick validation (valence checks, disconnected atoms)
   - Deep validation via backend RDKit API
   - Returns structured validation results

### UI Components

1. **Editor2D** (`components/Editor2D.tsx`)
   - 2D molecule editor using Konva
   - Grid snapping, hitbox detection
   - Atom placement, bond drawing, selection
   - Real-time updates from MoleculeStateEngine

2. **Viewer3D** (`components/Viewer3D.tsx`)
   - 3D molecule viewer using 3Dmol.js
   - Auto-syncs with 2D editor
   - Multiple rendering styles (stick, ballstick, sphere)

3. **ToolsSidebar** (`components/ToolsSidebar.tsx`)
   - Left sidebar with grouped tools
   - Basic Editing, Chemistry Tools, AI/Modeling sections
   - Undo/redo controls

4. **BottomDock** (`components/BottomDock.tsx`)
   - Context-aware bottom panel
   - Precision editor (atom properties)
   - Validation results
   - ML predictions

## Usage

### Basic Setup

```tsx
import Editor2D from './lab/components/Editor2D'
import Viewer3D from './lab/components/Viewer3D'
import ToolsSidebar from './lab/components/ToolsSidebar'
import BottomDock from './lab/components/BottomDock'
import { moleculeEngine } from './lab/engines/MoleculeStateEngine'
import { toolController, ToolMode } from './lab/engines/ToolController'
import { commandStack } from './lab/engines/CommandStack'
import { ValidationEngine } from './lab/engines/ValidationEngine'

function LabPage() {
  const [selectedElement, setSelectedElement] = useState<string | null>(null)
  const [selectedAtomId, setSelectedAtomId] = useState<string | null>(null)
  const [validationResult, setValidationResult] = useState(null)

  const handleValidate = async () => {
    const result = await ValidationEngine.validateCurrent()
    setValidationResult(result)
  }

  return (
    <div className="flex h-screen">
      <ToolsSidebar
        selectedElement={selectedElement}
        onElementSelect={setSelectedElement}
        onValidate={handleValidate}
      />
      <div className="flex-1 flex flex-col">
        <div className="flex-1 flex">
          <Editor2D
            selectedElement={selectedElement}
            onAtomSelect={setSelectedAtomId}
          />
          <Viewer3D />
        </div>
        <BottomDock
          selectedAtomId={selectedAtomId}
          validationResult={validationResult}
        />
      </div>
    </div>
  )
}
```

### Programmatic Usage

```tsx
import { moleculeEngine } from './lab/engines/MoleculeStateEngine'
import { commandStack, AddAtomCommand } from './lab/engines/CommandStack'

// Add an atom
const cmd = new AddAtomCommand('C', 100, 100)
commandStack.run(cmd)

// Get all atoms
const atoms = moleculeEngine.getAllAtoms()

// Serialize to SDF
const sdf = moleculeEngine.toSDF()

// Validate
const result = await ValidationEngine.validateCurrent()
```

## Dependencies

### Frontend

- **Konva** - For 2D canvas rendering
  ```bash
  npm install konva react-konva
  ```

- **3Dmol.js** - For 3D molecule visualization
  - Loaded via CDN (see Viewer3D.tsx)
  - Or install: `npm install 3dmol`

### Backend

- **RDKit** - Already installed via Conda
- **FastAPI** - Already installed

## Backend API

### Validation Endpoint

```bash
POST /api/validate
Content-Type: application/json

{
  "smiles": "CCO"
}

Response:
{
  "valid": true,
  "issues": [
    {
      "type": "warning",
      "message": "Atom has incomplete valence",
      "atomId": "123"
    }
  ]
}
```

### ML Prediction Endpoint

```bash
POST /api/ml/predict
Content-Type: application/json

{
  "atoms": [...],
  "bonds": [...]
}

Response:
{
  "logP": 0.22,
  "toxicity": 0.1,
  "stability": 0.87
}
```

## File Structure

```
frontend/src/lab/
├── engines/
│   ├── MoleculeStateEngine.ts    # Core molecule state
│   ├── ToolController.ts         # Tool mode management
│   ├── CommandStack.ts            # Undo/redo system
│   └── ValidationEngine.ts        # Validation logic
├── components/
│   ├── Editor2D.tsx              # 2D editor (Konva)
│   ├── Viewer3D.tsx              # 3D viewer (3Dmol.js)
│   ├── ToolsSidebar.tsx          # Left sidebar
│   └── BottomDock.tsx            # Bottom panel
└── README.md                      # This file

backend/
├── routes/
│   └── validate.py                # Validation API
└── services/
    └── property_predictor.py      # ML prediction service
```

## TODO / Next Steps

1. **Install Konva** for 2D editor:
   ```bash
   cd frontend
   npm install konva react-konva @types/konva
   ```

2. **Integrate with existing Lab page** - Replace or enhance current Lab.tsx

3. **Connect to backend** - Ensure validation and ML endpoints are working

4. **Add periodic table** - For element selection in ToolsSidebar

5. **Implement real SMILES/SDF conversion** - Currently using placeholders

6. **Load actual ML models** - Replace stubs in property_predictor.py

7. **Add unit tests** - For all engines

8. **Performance optimization** - For large molecules

## Notes

- All engines are singletons for easy access
- Commands are immutable - undo/redo is reliable
- Validation is two-tier: quick local + deep backend
- 2D and 3D viewers stay in sync via MoleculeStateEngine

