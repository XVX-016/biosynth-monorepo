"""
Property prediction model wrapper
"""
import os
import torch
from backend.ai.property_predictor import create_model
from backend.config import settings


class ModelLoader:
    """
    Loads and manages PropertyPredictor model
    """
    def __init__(self, weights_path: str = None):
        """
        Initialize model loader
        
        Args:
            weights_path: Path to model weights file (defaults to config)
        """
        if weights_path is None:
            weights_path = settings.MODEL_WEIGHTS_PATH
        self.weights_path = weights_path
        self.model = None
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self._load_model()
    
    def _load_model(self):
        """Load model from weights file"""
        if os.path.exists(self.weights_path):
            self.model = create_model()
            self.model.load_state_dict(torch.load(self.weights_path, map_location=self.device))
            self.model.eval()  # Set to evaluation mode
            self.model = self.model.to(self.device)
            print(f"Loaded model from {self.weights_path}")
        else:
            print(f"Warning: Model weights not found at {self.weights_path}")
            print("Using untrained model. Run train_property_predictor.py to train.")
            self.model = create_model()
            self.model.eval()
            self.model = self.model.to(self.device)
    
    def predict(self, features):
        """
        Predict molecular properties from features
        
        Args:
            features: Feature vector (list of floats, length 2048)
            
        Returns:
            dict with properties: stability, toxicity, solubility, bioavailability, novelty
        """
        if self.model is None:
            raise RuntimeError("Model not loaded")
        
        # Convert features to tensor
        if isinstance(features, list):
            x = torch.tensor(features, dtype=torch.float32).unsqueeze(0)  # Add batch dimension
        else:
            x = features
        
        x = x.to(self.device)
        
        # Run forward pass
        with torch.no_grad():
            y = self.model(x)
            y = y.cpu().numpy()[0]  # Remove batch dimension and convert to numpy
        
        return {
            "stability": float(y[0]),
            "toxicity": float(y[1]),
            "solubility": float(y[2]),
            "bioavailability": float(y[3]),
            "novelty": float(y[4]),
        }


class PropertyPredictor:
    """
    Wrapper for molecular property prediction model
    Maintains backward compatibility with existing code
    """
    def __init__(self, weights_path: str = None):
        if weights_path is None:
            weights_path = settings.MODEL_WEIGHTS_PATH
        self.loader = ModelLoader(weights_path)
    
    def predict(self, features):
        """
        Predict molecular properties from features
        
        Args:
            features: Feature vector (list or tensor)
            
        Returns:
            dict with properties: stability, toxicity, solubility, bioavailability, novelty
        """
        return self.loader.predict(features)


# Global predictor instance
_predictor = None
_last_predictions = None


def get_predictor():
    """Get or create global predictor instance"""
    global _predictor
    if _predictor is None:
        _predictor = PropertyPredictor()
    return _predictor


def get_last_predictions():
    """Get last predictions (for debugging)"""
    return _last_predictions


def set_last_predictions(predictions):
    """Set last predictions (for debugging)"""
    global _last_predictions
    _last_predictions = predictions

