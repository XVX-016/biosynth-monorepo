# Two-Table Architecture Implementation

Complete guide for the global + personal molecule library architecture.

## Overview

The system uses **two separate tables** for molecules:

1. **`public_molecules`** - Global, read-only library (curated molecules)
2. **`user_molecules`** - Personal user library (user-owned molecules)

## Database Schema

### public_molecules Table

```sql
CREATE TABLE public_molecules (
  id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
  name TEXT NOT NULL,
  formula TEXT,
  smiles TEXT,
  molfile TEXT,
  thumbnail_b64 TEXT,
  metadata JSONB DEFAULT '{}',
  created_at TIMESTAMPTZ DEFAULT now(),
  updated_at TIMESTAMPTZ DEFAULT now()
);
```

**RLS Policies:**
- ✅ **SELECT**: Anyone can read (public access)
- ❌ **INSERT/UPDATE/DELETE**: Only service role (admin operations)

### user_molecules Table

```sql
CREATE TABLE user_molecules (
  id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),
  user_id uuid REFERENCES auth.users NOT NULL,
  name TEXT NOT NULL,
  formula TEXT,
  smiles TEXT,
  molfile TEXT,
  thumbnail_b64 TEXT,
  metadata JSONB DEFAULT '{}',
  created_at TIMESTAMPTZ DEFAULT now(),
  updated_at TIMESTAMPTZ DEFAULT now()
);
```

**RLS Policies:**
- ✅ **SELECT**: Users can only see their own molecules
- ✅ **INSERT**: Users can create their own molecules
- ✅ **UPDATE**: Users can update their own molecules
- ✅ **DELETE**: Users can delete their own molecules

## Frontend Implementation

### Store Functions

#### Public Molecules (`frontend/src/lib/publicMoleculeStore.ts`)

```typescript
// List all public molecules
listPublicMolecules(): Promise<PublicMolecule[]>

// Search public molecules
searchPublicMolecules(query: string): Promise<PublicMolecule[]>

// Filter by element
filterPublicMoleculesByElement(element: string): Promise<PublicMolecule[]>

// Get by ID
getPublicMolecule(id: string): Promise<PublicMolecule | null>
```

#### User Molecules (`frontend/src/lib/userMoleculeStore.ts`)

```typescript
// Save user molecule
saveUserMolecule(userId: string, molecule: ...): Promise<string>

// List user molecules
listUserMolecules(userId: string): Promise<UserMolecule[]>

// Search user molecules
searchUserMolecules(userId: string, query: string): Promise<UserMolecule[]>

// Update user molecule
updateUserMolecule(userId: string, id: string, updates: ...): Promise<void>

// Delete user molecule
deleteUserMolecule(userId: string, id: string): Promise<void>

// Fork public molecule to user library
forkPublicMolecule(userId: string, publicId: string): Promise<string>
```

### Pages

#### Public Library (`/library/public`)

- Shows all public molecules
- Read-only access
- Search and filter support
- "Fork" button to copy to user library

#### User Library (`/library`)

- Shows user's personal molecules
- Full CRUD operations
- Search and filter support
- "Open in Lab" functionality

### Components

#### BarbellViewer (`frontend/src/components/BarbellViewer.tsx`)

- Three.js-based ball-and-stick 3D viewer
- Element-based coloring (CPK scheme)
- Interactive orbit controls
- Parses molfile format

#### MoleculeFilters (`frontend/src/components/MoleculeFilters.tsx`)

- Text search (name, formula, SMILES)
- Element filter
- Debounced input
- Clear button

## Backend Scripts

### Seed Public Molecules

```bash
# Generate SQL INSERT statements
python backend/scripts/seed_public_molecules.py
```

This generates SQL statements for 50 common molecules with:
- ✅ SMILES strings
- ✅ Formulas
- ✅ 3D molfiles (generated)
- ✅ 2D thumbnails (generated)

## Usage Patterns

### Viewing Public Molecules

```typescript
import { listPublicMolecules } from '../lib/publicMoleculeStore';

const molecules = await listPublicMolecules();
```

### Saving User Molecule

```typescript
import { saveUserMolecule } from '../lib/userMoleculeStore';

await saveUserMolecule(userId, {
  name: 'My Molecule',
  smiles: 'CCO',
  formula: 'C2H6O',
  molfile: '...',
  thumbnail_b64: '...',
});
```

### Forking Public Molecule

```typescript
import { forkPublicMolecule } from '../lib/userMoleculeStore';

// Creates a copy in user_molecules
const newId = await forkPublicMolecule(userId, publicMoleculeId);
```

## Benefits

### ✅ Clean Separation

- Public molecules: Curated, read-only
- User molecules: Personal, editable

### ✅ Security

- RLS policies enforce access control
- Users can't modify public molecules
- Users can only see their own molecules

### ✅ Scalability

- Public molecules: Shared across all users
- User molecules: Isolated per user
- Easy to add admin tools for public library

### ✅ Professional Pattern

- Same architecture as Notion, Figma, Google Drive
- Industry-standard approach
- Easy to understand and maintain

## Migration from Single Table

If you have existing molecules in a single `molecules` table:

1. **Backup existing data**
2. **Create new tables** (`public_molecules`, `user_molecules`)
3. **Migrate data**:
   - Molecules with `user_id` → `user_molecules`
   - Molecules without `user_id` → `public_molecules`
4. **Update frontend** to use new stores
5. **Test thoroughly**

## Next Steps

1. ✅ Tables created (you mentioned already done)
2. ✅ RLS policies configured
3. ✅ Frontend stores created
4. ✅ Components created (BarbellViewer, MoleculeFilters)
5. ✅ Pages created (PublicLibrary)
6. ✅ Lab page updated to save to user_molecules
7. ⏳ Seed public_molecules with 50 molecules
8. ⏳ Add routes to App.tsx
9. ⏳ Test fork functionality

## Testing Checklist

- [ ] Public molecules visible to all users
- [ ] User molecules only visible to owner
- [ ] Fork functionality works
- [ ] Search and filters work
- [ ] Barbell viewer renders correctly
- [ ] Lab saves to user_molecules
- [ ] RLS policies prevent unauthorized access

