# Phase 8: Conformer Generator

## Overview

Phase 8 provides interfaces and mock implementations for conformer generation.

## Structure

```
backend/src/conformers/
├── __init__.py              # Module exports
├── conformer_generator.py   # Mock conformer generator
└── etkdg.py                # ETKDG generator (placeholder)
```

## Features

### Conformer Generator (`conformer_generator.py`)
- `generate_conformers()` - Generate n conformers
- Returns sorted list (by energy)
- Mock 3D coordinates

### ETKDG Generator (`etkdg.py`)
- Placeholder for RDKit ETKDG integration
- Falls back to mock generator

## Usage

```python
from src.conformers import ConformerGenerator

generator = ConformerGenerator()
conformers = generator.generate_conformers("CCO", n=10)
print(f"Generated {len(conformers)} conformers")
```

## API Endpoints

- `POST /api/conformers/generate` - Generate conformers

## Future Integration

To integrate real conformer generation:
1. Use RDKit's ETKDG implementation
2. Add conformer optimization
3. Replace mock coordinates with real 3D structures

