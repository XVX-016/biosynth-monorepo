# BioSynth AI Monorepo

Full-stack molecular design system with React + R3F frontend, TypeScript engine, and FastAPI + PyTorch backend.

## Architecture

```
biosynth-monorepo/
├── frontend/          # React + TypeScript + Tailwind + Framer Motion + R3F
├── backend/           # FastAPI + PyTorch + RDKit
└── packages/
    └── engine/        # Pure TypeScript molecular engine
```

## Tech Stack

### Frontend
- **React 18.2.0** + **TypeScript 5.3.3**
- **Vite 5.0.8** (build tool)
- **TailwindCSS 3.4.1** (aluminium design system)
- **Framer Motion 11.0.3** (animations)
- **react-three-fiber 8.15.9** + **drei 9.101.3** (3D rendering)
- **Zustand 4.5.0** (state management)
- **React Router 6.22.0**

### Engine (packages/engine)
- **TypeScript 5.3.3** (pure, no React/Node dependencies)
- **Vitest** (testing)

### Backend
- **Python 3.10**
- **FastAPI 0.110.0**
- **PyTorch 2.2.0**
- **RDKit 2023.09.2**
- **ONNX Runtime 1.17.0**
- **Pytest** (testing)

## Getting Started

### Prerequisites

- Node.js 18+
- Python 3.10+
- npm or yarn

### Installation

1. **Bootstrap the monorepo:**
   ```bash
   npm install
   ```

2. **Set up backend Python environment:**
   ```bash
   cd backend
   python -m venv .venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate
   pip install -r requirements.txt
   cd ..
   ```

### Development

**Run all services concurrently:**
```bash
npm run dev
```

**Run individually:**
```bash
# Frontend (http://localhost:5173)
npm run dev:frontend

# Backend (http://localhost:8000)
npm run dev:backend
```

### Testing

**Run all tests:**
```bash
./run_all_tests.sh
```

**Run tests individually:**
```bash
# Engine tests
cd packages/engine && npm test

# Frontend tests
cd frontend && npm test

# Backend tests
cd backend && pytest
```

### Building

**Build all packages:**
```bash
npm run build
```

## Project Structure

### Frontend (`frontend/`)
- `src/pages/` - Page components
- `src/components/` - Reusable UI components
- `src/components/r3f/` - React Three Fiber 3D components
- `src/store/` - Zustand state management
- Uses aluminium design palette from `tailwind.config.js`

### Engine (`packages/engine/`)
- `src/types.ts` - Core types (Atom, Bond, Element)
- `src/MoleculeGraph.ts` - Molecular graph data structure
- `src/LayoutEngine.ts` - 3D positioning and force field optimization
- `src/UndoStack.ts` - Undo/redo functionality
- Pure TypeScript, no external dependencies (except testing)

### Backend (`backend/`)
- `app.py` - FastAPI application
- `models/` - ML models (to be implemented)
- `utils/` - RDKit utilities (to be implemented)
- `weights/` - Model weights (gitignored)

## Design System

The UI follows the **BioSynth Aluminium** design system:

- **Aluminium colors**: `#F5F5F7` (light), `#E5E7EA` (default), `#C9CCD1` (dark)
- **Accent colors**: `#4EA7FF` (blue), `#6EE787` (green), `#FF6B6B` (red)
- **Text colors**: `#1A1A1C` (primary), `#5C5C5F` (secondary), `#9A9A9E` (tertiary)
- **Border radius**: `18px` for large elements
- **Shadows**: `elev-1` for elevation

## API Endpoints

### Backend (`http://localhost:8000`)

- `GET /` - API information
- `GET /health` - Health check
- `POST /predict` - Predict molecular properties from SMILES
  ```json
  {
    "smiles": "CCO"
  }
  ```

## Development Guidelines

See `cursor.json` for detailed coding standards:

- **Frontend**: React + TypeScript + Tailwind + Framer Motion, R3F for 3D
- **Engine**: Pure TypeScript, deterministic functions, unit-testable
- **Backend**: FastAPI, modular ML models, RDKit in utils
- **Testing**: Vitest (frontend/engine), Pytest (backend)
- **Type Safety**: Strong typing required, no `any` types

## Next Steps

1. Implement full molecular engine features (force field, SMILES export)
2. Add ML model training scripts
3. Implement molecule generator (Transformer-based)
4. Add ONNX export for production inference
5. Set up CI/CD pipeline
6. Add comprehensive documentation

## License

Private - All Rights Reserved

