# Next Steps - Wiring LabV2 Buttons

## Quick Guide to Complete the Interactive Loop

### 1. Wire the Optimize Button

**File:** `frontend/src/components/LabV2/FloatingToolbar.tsx`

**Current:**
```tsx
<button 
  className="lab-btn lab-btn-primary" 
  onClick={() => dispatch({ type: "SET_BUSY", payload: true })}
  disabled={state.busy}
>
  Optimize
</button>
```

**Replace with:**
```tsx
<button 
  className="lab-btn lab-btn-primary" 
  onClick={async () => {
    dispatch({ type: "SET_BUSY", payload: true });
    try {
      const { MLAPI } = await import("../../api/ml");
      const payload = {
        atoms: state.atoms.map(a => ({
          id: a.id,
          element: a.element,
          x: a.position[0],
          y: a.position[1],
          z: a.position[2] || 0
        })),
        bonds: state.bonds.map(b => ({
          a: b.a,
          b: b.b,
          order: b.order || 1
        }))
      };
      const res = await MLAPI.optimize(payload);
      dispatch({ type: "APPLY_ML", payload: { correctedAtoms: res.correctedAtoms } });
    } catch (e) {
      console.error("Optimize failed", e);
    } finally {
      dispatch({ type: "SET_BUSY", payload: false });
    }
  }}
  disabled={state.busy}
>
  Optimize
</button>
```

---

### 2. Wire the Predict Button

**File:** `frontend/src/components/LabV2/FloatingToolbar.tsx`

**Current:**
```tsx
<button 
  className="lab-btn" 
  onClick={() => dispatch({ type: "SET_BUSY", payload: true })}
  disabled={state.busy}
>
  Predict
</button>
```

**Replace with:**
```tsx
<button 
  className="lab-btn" 
  onClick={async () => {
    dispatch({ type: "SET_BUSY", payload: true });
    try {
      const response = await fetch("http://localhost:8000/lab/predict-bonds", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          atoms: state.atoms.map(a => ({
            id: a.id,
            element: a.element,
            position: a.position
          }))
        })
      });
      const data = await response.json();
      
      // Add predicted bonds to state
      const newBonds = data.bonds.map(b => ({
        id: b.id,
        a: b.a,
        b: b.b,
        order: b.order
      }));
      
      dispatch({ 
        type: "LOAD_MOLECULE", 
        payload: { 
          atoms: state.atoms, 
          bonds: newBonds 
        } 
      });
    } catch (e) {
      console.error("Predict failed", e);
    } finally {
      dispatch({ type: "SET_BUSY", payload: false });
    }
  }}
  disabled={state.busy}
>
  Predict
</button>
```

---

### 3. Add Export Button

**File:** `frontend/src/components/LabV2/FloatingToolbar.tsx`

**Add after the Clear button:**
```tsx
<button 
  className="lab-btn" 
  onClick={async () => {
    try {
      const response = await fetch("http://localhost:8000/lab/export-molecule", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          atoms: state.atoms.map(a => ({
            id: a.id,
            element: a.element,
            x: a.position[0],
            y: a.position[1],
            z: a.position[2] || 0
          })),
          bonds: state.bonds.map(b => ({
            a: b.a,
            b: b.b,
            order: b.order || 1
          })),
          format: "molforge"  // or "mol" or "pdb"
        })
      });
      const data = await response.json();
      
      // Download file
      const blob = new Blob([data.data], { type: "application/json" });
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = `molecule.${data.format}`;
      a.click();
      URL.revokeObjectURL(url);
    } catch (e) {
      console.error("Export failed", e);
    }
  }}
>
  Export
</button>
```

---

### 4. Test the Complete Flow

1. **Start backend:**
   ```bash
   cd backend
   python app.py
   ```

2. **Start frontend:**
   ```bash
   cd frontend
   npm run dev
   ```

3. **Navigate to:** `http://localhost:5173/labv2`

4. **Test workflow:**
   - Click "Add Atom"
   - Click canvas to place 2-3 carbon atoms
   - Click "Predict" → bonds should appear
   - Click "Optimize" → atoms should adjust positions
   - Click "Export" → download .molforge file

---

### 5. Optional: Add Keyboard Shortcuts

**File:** `frontend/src/components/LabV2/LabV2Page.tsx`

**Add inside the main div:**
```tsx
useEffect(() => {
  const handleKeyDown = (e: KeyboardEvent) => {
    if (e.ctrlKey && e.key === 'z') {
      e.preventDefault();
      undo();
    }
    if (e.ctrlKey && e.key === 'y') {
      e.preventDefault();
      redo();
    }
    if (e.key === 'Delete' || e.key === 'Backspace') {
      if (state.selection) {
        dispatch({ type: "DELETE_ATOM", payload: state.selection });
      }
    }
  };
  
  window.addEventListener('keydown', handleKeyDown);
  return () => window.removeEventListener('keydown', handleKeyDown);
}, [undo, redo, state.selection, dispatch]);
```

---

### 6. Verify Backend is Running

**Test predict-bonds endpoint:**
```bash
curl -X POST http://localhost:8000/lab/predict-bonds \
  -H "Content-Type: application/json" \
  -d '{
    "atoms": [
      {"id": "0", "element": "C", "position": [0, 0, 0]},
      {"id": "1", "element": "H", "position": [1, 0, 0]}
    ]
  }'
```

**Expected response:**
```json
{
  "bonds": [
    {"id": "b_0_1", "a": "0", "b": "1", "order": 1}
  ]
}
```

---

## Troubleshooting

### CORS Errors

If you see CORS errors, ensure backend has CORS middleware:

**File:** `backend/app.py`

```python
from fastapi.middleware.cors import CORSMiddleware

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
```

### Import Errors

If `MLAPI` import fails, check that `frontend/src/api/ml.ts` exists and exports `MLAPI`.

### Type Errors

If TypeScript complains about atom properties, ensure `EditorContext` Atom type matches:

```tsx
export type Atom = { 
  id: string; 
  element: string; 
  position: [number, number, number] 
};
```

---

## Quick Win Checklist

- [ ] Wire Optimize button
- [ ] Wire Predict button  
- [ ] Add Export button
- [ ] Test full workflow
- [ ] Add keyboard shortcuts (optional)
- [ ] Verify CORS is enabled

---

**Estimated time:** 15-20 minutes

**Result:** Fully functional molecular editor with ML-powered bond prediction!
