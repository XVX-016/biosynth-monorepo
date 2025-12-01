# API Fixes Applied

## Issues Fixed

### 1. PyTorch Import Error (Line 70)
**Problem**: Server was catching OSError during PyTorch import, setting TORCH_AVAILABLE=False
**Fix**: Added retry logic in `featurize_smiles()` to re-import PyTorch if initial import failed
**Location**: `backend/ml/featurize.py`

### 2. Batch Endpoint 422 Error (Line 72)
**Problem**: Endpoint expected `molecules` but test script sent `inputs`
**Fix**: Updated `BatchPredictRequest` to accept both `molecules` and `inputs` formats
**Location**: `backend/api/predict.py`

### 3. Attention Map Processing
**Problem**: Potential None/tensor issues in attention processing
**Fix**: Already fixed in previous commit with TORCH_AVAILABLE checks

## Testing

After restarting the server, test with:
```bash
python backend/scripts/test_api_endpoints.py
```

## Server Restart Required

The server needs to be restarted to pick up the PyTorch import fixes:
```bash
cd backend
$env:PYTHONPATH = "C:\Computing\MolForge"
python -m uvicorn app:app --host 0.0.0.0 --port 8000
```

