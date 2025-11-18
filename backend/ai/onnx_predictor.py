"""
ONNX-based predictor for accelerated inference
"""
import os
import numpy as np
from backend.config import settings

try:
    import onnxruntime as ort
    ONNXRUNTIME_AVAILABLE = True
except ImportError:
    ONNXRUNTIME_AVAILABLE = False
    ort = None


class ONNXPredictor:
    """
    ONNX-based predictor for fast inference
    """
    def __init__(self, onnx_path: str = None):
        """
        Initialize ONNX predictor
        
        Args:
            onnx_path: Path to ONNX model file (defaults to config)
        """
        if onnx_path is None:
            onnx_path = settings.ONNX_MODEL_PATH
        if not ONNXRUNTIME_AVAILABLE:
            raise ImportError(
                "onnxruntime is not installed. Install it with: pip install onnxruntime"
            )
        
        self.onnx_path = onnx_path
        
        if not os.path.exists(onnx_path):
            raise FileNotFoundError(
                f"ONNX model not found at {onnx_path}. "
                "Run export_to_onnx.py to create ONNX model."
            )
        
        # Create inference session
        self.session = ort.InferenceSession(onnx_path)
        self.input_name = self.session.get_inputs()[0].name
        self.output_name = self.session.get_outputs()[0].name
        
        print(f"Loaded ONNX model from {onnx_path}")
    
    def predict(self, features):
        """
        Predict molecular properties from features using ONNX
        
        Args:
            features: Feature vector (list of floats, length 2048)
            
        Returns:
            dict with properties: stability, toxicity, solubility, bioavailability, novelty
        """
        # Convert features to numpy array
        if isinstance(features, list):
            x = np.array(features, dtype=np.float32).reshape(1, -1)  # Add batch dimension
        else:
            x = np.array(features, dtype=np.float32)
            if x.ndim == 1:
                x = x.reshape(1, -1)
        
        # Run inference
        result = self.session.run([self.output_name], {self.input_name: x})
        y = result[0][0]  # Remove batch dimension
        
        return {
            "stability": float(y[0]),
            "toxicity": float(y[1]),
            "solubility": float(y[2]),
            "bioavailability": float(y[3]),
            "novelty": float(y[4]),
        }


# Global ONNX predictor instance
_onnx_predictor = None


def get_onnx_predictor(onnx_path: str = None):
    """
    Get or create global ONNX predictor instance
    
    Args:
        onnx_path: Path to ONNX model file (defaults to config)
    
    Returns:
        ONNXPredictor instance if onnxruntime is available, None otherwise
    """
    if onnx_path is None:
        onnx_path = settings.ONNX_MODEL_PATH
    if not ONNXRUNTIME_AVAILABLE:
        return None
    
    global _onnx_predictor
    if _onnx_predictor is None:
        _onnx_predictor = ONNXPredictor(onnx_path)
    return _onnx_predictor

