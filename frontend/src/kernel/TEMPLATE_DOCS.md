# Template System Overview

## Files
- `templates/*.json` — canonical molecule templates
- `templateLoader.ts` — loading & placing templates
- `TemplatePicker.tsx` — UI for selecting templates

## Usage

### Load + place:
```ts
const t = loadTemplate("water");
placeTemplate(t, {x: 0, y: 0, z: 0});
```

### UI Component:
```tsx
<TemplatePicker />
```

## Graph Integration

Template placement directly modifies the store's MoleculeGraph:

- Atoms converted 1:1
- Bonds created from local indices
- Returns `createdAtomIds[]` for selection logic
- Automatically pushes to history store

## Template Format

Templates are JSON files with:
- `atoms`: Array of `{ element: string, coords: { x, y, z } }`
- `bonds`: Array of `{ a: number, b: number, order: number }` (indices into atoms array)

