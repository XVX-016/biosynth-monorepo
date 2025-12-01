# Python 3.10 Setup Guide

## ✅ PyTorch Installation Complete

This project now uses **Python 3.10** with a clean virtual environment.

### Environment Setup

```powershell
cd backend
py -3.10 -m venv venv
.\venv\Scripts\activate
```

### Install PyTorch (CPU)

```powershell
pip install torch==2.2.2 --index-url https://download.pytorch.org/whl/cpu
pip install torch-geometric==2.5.2
pip install torch-scatter torch-sparse torch-cluster torch-spline-conv -f https://data.pyg.org/whl/torch-2.2.0+cpu.html
```

### Install Dependencies

```powershell
pip install -r requirements.txt
pip install "numpy<2.0"  # Fix NumPy 2.x compatibility
```

### Verify Installation

```powershell
python -c "import torch; import torch_geometric; print('PyTorch:', torch.__version__); print('PyG:', torch_geometric.__version__)"
```

### Start Server

```powershell
$env:PYTHONPATH = "C:\Computing\MolForge"
python -m uvicorn app:app --host 0.0.0.0 --port 8000
```

## Fixed Issues

1. ✅ **PyTorch DLL Error**: Fixed by using Python 3.10 + clean venv
2. ✅ **Batch Endpoint 422**: Fixed by handling both `molecules` and `inputs` formats
3. ✅ **NumPy Compatibility**: Downgraded to NumPy < 2.0 for PyTorch 2.2.2

## Notes

- Python 3.10 is required (not 3.11 or 3.14)
- Use the venv in `backend/venv/`
- Always activate venv before running server: `.\venv\Scripts\activate`

