# backend/ai/__init__.py
from .featurizer import featurize_smiles, validate_smiles
# Lazy import predictor to avoid torch DLL issues on startup
try:
    from .predictor import get_predictor, PropertyPredictor
except (ImportError, RuntimeError):
    # If torch is not available, provide fallback
    def get_predictor():
        raise RuntimeError("PyTorch not available. Install with: pip install torch --index-url https://download.pytorch.org/whl/cpu")
    PropertyPredictor = None
from .onnx_predictor import get_onnx_predictor

__all__ = [
    "featurize_smiles",
    "validate_smiles",
    "get_predictor",
    "PropertyPredictor",
    "get_onnx_predictor",
]

