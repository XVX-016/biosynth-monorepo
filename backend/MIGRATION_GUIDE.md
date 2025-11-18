# Python 3.14 to 3.11 Migration Guide

This guide documents the migration from Python 3.14 to Python 3.11 and the fixes applied to resolve installation issues.

## Changes Made

### 1. Updated Package Versions (`requirements.txt`)
- **PyTorch**: Changed from `2.9.1` (doesn't exist) to `2.4.1` (stable for Python 3.11)
- **NumPy**: Changed from `2.1.2` to `1.26.4` (stable for Python 3.11)
- **Pydantic**: Changed from `2.12.4` to `2.9.2` (compatible with other dependencies)
- **Transformers**: Changed from `4.50.0` to `4.41.2` (stable version)
- **Tokenizers**: Changed from `0.21.4` to `0.19.1` (matching transformers)
- **Pillow**: Changed from `11.3.0` to `10.4.0` (stable version)
- **RDKit**: Uncommented `rdkit-pypi==2024.03.4` (now available for Python 3.11)
- **ONNXRuntime**: Uncommented `onnxruntime==1.19.2` (now available for Python 3.11)
- **Added**: `pandas==2.2.2` and `tqdm==4.66.4` (required by training scripts)

### 2. Updated Dockerfile
- Changed base image from `python:3.10-slim` to `python:3.11-slim`

### 3. Code Refactoring
- **`onnx_predictor.py`**: Added graceful handling for missing `onnxruntime` with optional import
- **`app.py`**: 
  - Added environment variable loading with `python-dotenv`
  - Improved error handling for ONNX predictor fallback

### 4. Environment Variables
- Created `env.example` file with all required environment variables
- Added `load_dotenv()` to `app.py` for automatic `.env` file loading

### 5. Documentation
- Created `STRUCTURE.md` documenting the backend architecture
- Created this migration guide

## Setup Instructions

### Step 1: Create Fresh Virtual Environment

**Windows:**
```powershell
# Remove old virtual environment if it exists
Remove-Item -Recurse -Force .venv -ErrorAction SilentlyContinue

# Create new virtual environment with Python 3.11
python3.11 -m venv .venv

# Activate virtual environment
.venv\Scripts\Activate.ps1
```

**Linux/Mac:**
```bash
# Remove old virtual environment if it exists
rm -rf .venv

# Create new virtual environment with Python 3.11
python3.11 -m venv .venv

# Activate virtual environment
source .venv/bin/activate
```

### Step 2: Verify Python Version
```bash
python --version
# Should output: Python 3.11.x
```

### Step 3: Upgrade pip
```bash
python -m pip install --upgrade pip
```

### Step 4: Install Dependencies
```bash
cd backend
pip install -r requirements.txt
```

**Note**: If you encounter issues with specific packages:
- RDKit installation can be slow - be patient
- PyTorch installation may require specific CUDA versions if using GPU
- If onnxruntime fails, you can comment it out (the code handles it gracefully)

### Step 5: Set Up Environment Variables
```bash
# Copy the example file
cp env.example .env

# Edit .env with your configuration
# At minimum, set:
# - DATABASE_URL (default SQLite is fine for development)
# - REDIS_URL (if using background jobs)
```

### Step 6: Verify Installation
```bash
# Test imports
python -c "import torch; import rdkit; import onnxruntime; print('All imports successful!')"

# Run the application
uvicorn app:app --reload
```

## Troubleshooting

### Issue: RDKit installation fails
**Solution**: RDKit can be tricky. Try:
```bash
pip install rdkit-pypi==2024.03.4 --no-cache-dir
```

### Issue: PyTorch installation fails
**Solution**: Install PyTorch from the official site:
```bash
# CPU only
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu

# CUDA (adjust version as needed)
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121
```

### Issue: ONNXRuntime installation fails
**Solution**: The code handles missing onnxruntime gracefully. You can:
1. Skip it if not needed (ONNX features will fallback to PyTorch)
2. Install CPU version: `pip install onnxruntime`
3. Install GPU version: `pip install onnxruntime-gpu`

### Issue: Import errors after installation
**Solution**: 
1. Verify you're in the correct virtual environment
2. Check Python version: `python --version`
3. Reinstall packages: `pip install --force-reinstall -r requirements.txt`

## Testing the Migration

After setup, test the key components:

```bash
# Test database connection
python -c "from backend.db import init_db; init_db(); print('Database OK')"

# Test RDKit featurization
python -c "from backend.ai.featurizer import featurize_smiles; print(featurize_smiles('CCO'))"

# Test model loading
python -c "from backend.ai.property_predictor import create_model; m = create_model(); print('Model OK')"

# Run tests
pytest
```

## Next Steps

1. **Train a model** (if needed):
   ```bash
   python -m backend.models.train_property_predictor
   ```

2. **Export to ONNX** (optional, for faster inference):
   ```bash
   python -m backend.models.export_to_onnx
   ```

3. **Start the server**:
   ```bash
   uvicorn app:app --reload --host 0.0.0.0 --port 8000
   ```

4. **Start background worker** (if using Redis):
   ```bash
   python -m backend.jobs.worker
   ```

## Rollback Plan

If you need to rollback:
1. Restore the old `requirements.txt` from git
2. Use Python 3.14 virtual environment
3. Note: Some packages may not be available for Python 3.14

