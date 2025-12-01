# Integration Tests

## Overview

Integration tests verify that Phases 7-9 work correctly together and with previous phases.

## Test Files

- `test_phase7_integration.py` - Search and screening engine tests
- `test_phase8_integration.py` - Conformer generator tests
- `test_phase9_integration.py` - Multi-agent orchestrator tests
- `test_attention_integration.py` - Phase 5 attention mapping tests (existing)

## Running Tests

```bash
# Run all integration tests
pytest backend/tests/test_phase*_integration.py -v

# Run specific phase
pytest backend/tests/test_phase7_integration.py -v
pytest backend/tests/test_phase8_integration.py -v
pytest backend/tests/test_phase9_integration.py -v

# Run with coverage
pytest backend/tests/test_phase*_integration.py --cov=backend/src --cov-report=html
```

## Test Coverage

### Phase 7 Tests
- ✅ Similarity search
- ✅ Substructure search
- ✅ Library loader
- ✅ Screening pipeline
- ✅ API route imports

### Phase 8 Tests
- ✅ Conformer generation
- ✅ ETKDG fallback
- ✅ Invalid SMILES handling
- ✅ API route imports

### Phase 9 Tests
- ✅ Agent registration
- ✅ Task submission
- ✅ Task execution (async)
- ✅ Workflow execution
- ✅ Task status tracking
- ✅ All agent types (predictor, screening, qm, md)
- ✅ API route imports

## Notes

- Some tests may skip if dependencies are not available
- Mock implementations are tested (no real QM/MD required)
- Async tests use `pytest.mark.asyncio`

