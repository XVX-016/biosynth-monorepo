# Phase 8: Molecular Dynamics Engine

## Overview

Phase 8 provides interfaces and mock implementations for MD simulations. All results are fake but follow realistic patterns.

## Structure

```
backend/src/md/
├── __init__.py              # Module exports
├── md_interfaces.py         # Protocols for MD engines
├── md_engine.py            # Mock MD engine
├── forcefields.py          # Force field implementations
└── integrators.py          # MD integrators
```

## Features

### MD Engine (`md_engine.py`)
- `simulate()` - Mock MD simulation
- Returns trajectory with frames
- Computes energy, temperature, etc.

### Force Fields (`forcefields.py`)
- `SimpleForceField` - Mock harmonic + LJ-like force field
- `ForceField` protocol for real implementations

### Integrators (`integrators.py`)
- `VerletIntegrator` - Mock Velocity Verlet integrator
- `Integrator` protocol for real implementations

## Usage

```python
from src.md import MDEngine

engine = MDEngine()
result = await engine.simulate(
    smiles="CCO",
    steps=1000,
    timestep=0.001,
    temperature=300.0
)
print(f"Trajectory frames: {len(result.trajectory)}")
```

## API Endpoints

- `POST /api/md/simulate` - Run MD simulation

## Future Integration

To integrate real MD engines:
1. Implement `MDEngineProtocol` in engine classes
2. Add actual MD package dependencies (OpenMM, GROMACS, etc.)
3. Replace mock methods with real simulations

