# History System

## Phase 4 Complete ✅

Undo/Redo system for molecule editor operations.

## Structure

```
history/
├── Command.ts          # Command interface and implementations
├── HistoryManager.ts   # History stack management
└── index.ts           # Public API
```

## Commands

All molecule operations are commands:
- `AddAtomCommand` - Add an atom
- `RemoveAtomCommand` - Remove an atom (and its bonds)
- `AddBondCommand` - Add a bond
- `RemoveBondCommand` - Remove a bond
- `MoveAtomCommand` - Move an atom
- `UpdateAtomCommand` - Update atom properties
- `UpdateBondCommand` - Update bond properties
- `ClearMoleculeCommand` - Clear entire molecule

## Usage

```typescript
import { HistoryManager, AddAtomCommand } from '@/lib/molecule/history'

const history = new HistoryManager()
const molecule = new Molecule()

// Execute command
const command = new AddAtomCommand('atom1', 'C', [0, 0, 0])
history.execute(molecule, command)

// Undo
history.undo(molecule)

// Redo
history.redo(molecule)
```

## Features

- ✅ Command pattern for all operations
- ✅ Full state restoration on undo
- ✅ History size limit (100 commands)
- ✅ Keyboard shortcuts (Ctrl+Z, Ctrl+Shift+Z)
- ✅ UI state tracking (canUndo, canRedo)

