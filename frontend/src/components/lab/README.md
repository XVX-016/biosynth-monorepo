# Lab Foundation Architecture

This directory contains the new foundation architecture for the Lab page, built according to the specification for a clean, maintainable 3D molecule editor.

## Architecture Overview

### Core Components

1. **Types** (`../../types/molecule.ts`)
   - Canonical `Atom`, `Bond`, and `Molecule` types
   - Simple array-based structure (not Maps)
   - `ToolName` enum for tool types

2. **Store** (`../../store/labStore.ts`)
   - Zustand store as single source of truth
   - Undo/redo support built-in
   - All molecule mutations go through store actions
   - DevTools integration for debugging

3. **3D Rendering** (`SceneRoot.tsx`, `MoleculeScene.tsx`, `AtomMesh.tsx`, `BondMesh.tsx`)
   - Clean separation: rendering is read-only from store
   - Non-instanced meshes for now (can upgrade to InstancedMesh later)
   - Metallic/ivory theme matching existing design

4. **Tools** (`../../tools/`)
   - Clean tool interface
   - Isolated tool implementations (select, add_atom, bond, delete)
   - Tools dispatch actions to store, don't mutate directly

5. **Bonding Engine** (`../../utils/bondingEngine.ts`)
   - Valence-based bonding rules
   - Distance threshold checks
   - Prevents overbonding

6. **Serialization** (`../../utils/serialize.ts`)
   - JSON serialization/deserialization
   - Ready for MOL/SDF export later

## Usage

### Basic Example

```tsx
import LabViewer from '../components/lab/LabViewer'
import LabToolPanel from '../components/lab/LabToolPanel'
import { useLabStore } from '../store/labStore'

function MyLabPage() {
  return (
    <div>
      <LabToolPanel />
      <LabViewer />
    </div>
  )
}
```

### Programmatic Usage

```tsx
import { useLabStore } from '../store/labStore'

// Add an atom
const store = useLabStore.getState()
store.addAtom('C', [0, 0, 0])

// Add a bond
store.addBond(atomId1, atomId2, 1)

// Undo/redo
store.undo()
store.redo()

// Load a molecule
store.loadMolecule({
  atoms: [...],
  bonds: [...],
  name: 'My Molecule'
})
```

## Integration with Existing Code

The new foundation runs alongside the existing `MoleculeGraph`-based system. To migrate:

1. Create adapter functions to convert between `MoleculeGraph` and `Molecule` types
2. Gradually replace `useMoleculeStore` with `useLabStore` in components
3. Update `MoleculeViewer` to use `SceneRoot` and `MoleculeScene`

## Next Steps

- [ ] Convert to InstancedMesh for performance
- [ ] Add move/drag tool implementation
- [ ] Integrate with backend API for saving
- [ ] Add molecule templates loader
- [ ] Implement Explore Mode (physics simulation)
- [ ] Add molecule scoring/validation

## Testing

Test the foundation:

1. Navigate to `/lab-new` (if route added)
2. Click "Add Atom" tool
3. Select an element (C, H, O, etc.)
4. Click on canvas to place atoms
5. Use "Bond" tool to connect atoms
6. Test undo/redo functionality

## File Structure

```
frontend/src/
├── types/
│   └── molecule.ts          # Canonical types
├── store/
│   └── labStore.ts          # Zustand store
├── components/lab/
│   ├── SceneRoot.tsx        # Canvas wrapper
│   ├── MoleculeScene.tsx    # Scene container
│   ├── AtomMesh.tsx          # Atom rendering
│   ├── BondMesh.tsx          # Bond rendering
│   ├── ToolHandler.tsx       # Tool event routing
│   ├── LabToolPanel.tsx      # UI controls
│   └── LabViewer.tsx         # Main viewer component
├── tools/
│   ├── toolInterface.ts      # Tool interface
│   ├── selectTool.ts
│   ├── addAtomTool.ts
│   ├── bondTool.ts
│   ├── deleteTool.ts
│   └── index.ts
├── utils/
│   ├── bondingEngine.ts     # Auto-bonding logic
│   └── serialize.ts          # Serialization
└── data/
    └── elementPalette.ts     # Element colors/radii
```

