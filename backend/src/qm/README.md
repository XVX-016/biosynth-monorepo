# Phase 8: Quantum Chemistry Engine

## Overview

Phase 8 provides interfaces and mock implementations for quantum chemistry calculations. All results are fake but follow realistic patterns.

## Structure

```
backend/src/qm/
├── __init__.py              # Module exports
├── qm_interfaces.py         # Protocols for QM engines
├── qm_engine.py            # Mock QM engine
├── psi4_wrapper.py         # Psi4 interface (placeholder)
└── xtb_wrapper.py          # xTB interface (placeholder)
```

## Features

### QM Engine (`qm_engine.py`)
- `compute_energy()` - Mock electronic energy calculation
- `optimize_geometry()` - Mock geometry optimization
- `compute_properties()` - Mock property calculation (dipole, charges)

### Interfaces (`qm_interfaces.py`)
- `QMEngineProtocol` - Protocol for real QM engines
- `QMResult` - Result dataclass with energy, coordinates, forces, etc.

### Wrappers
- `Psi4Wrapper` - Placeholder for Psi4 integration
- `XTBWrapper` - Placeholder for xTB integration

## Usage

```python
from src.qm import QMEngine

engine = QMEngine()
result = await engine.compute_energy("CCO", method="B3LYP", basis="6-31G")
print(f"Energy: {result.energy} Hartree")
```

## API Endpoints

- `POST /api/qm/energy` - Compute QM energy
- `POST /api/qm/optimize` - Optimize geometry

## Future Integration

To integrate real QM engines:
1. Implement `QMEngineProtocol` in wrapper classes
2. Add actual QM package dependencies
3. Replace mock methods with real calculations

