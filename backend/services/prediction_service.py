# backend/services/prediction_service.py
"""
Service layer for molecular property prediction
"""
from typing import Dict, Optional
from backend.ai.featurizer import featurize_smiles, validate_smiles
from backend.ai.predictor import get_predictor
from backend.ai.onnx_predictor import get_onnx_predictor
from backend.core.exceptions import InvalidSMILESError, PredictionError


class PredictionService:
    """Service for molecular property prediction"""
    
    @staticmethod
    def predict_properties(smiles: str, use_onnx: bool = False) -> Dict[str, float]:
        """
        Predict molecular properties from SMILES string
        
        Args:
            smiles: SMILES string
            use_onnx: Whether to use ONNX model (faster)
            
        Returns:
            Dictionary of predicted properties
            
        Raises:
            InvalidSMILESError: If SMILES is invalid
            PredictionError: If prediction fails
        """
        # Validate SMILES
        if not validate_smiles(smiles):
            raise InvalidSMILESError(f"Invalid SMILES string: {smiles}")
        
        # Featurize
        features = featurize_smiles(smiles)
        if features is None:
            raise PredictionError(f"Failed to featurize SMILES: {smiles}")
        
        # Predict
        try:
            if use_onnx:
                predictor = get_onnx_predictor()
                if predictor is None:
                    # Fallback to regular predictor
                    predictor = get_predictor()
            else:
                predictor = get_predictor()
            
            properties = predictor.predict(features)
            return properties
        except Exception as e:
            raise PredictionError(f"Prediction failed: {str(e)}")

