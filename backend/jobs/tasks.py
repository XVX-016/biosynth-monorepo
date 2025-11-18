# backend/jobs/tasks.py
"""
Background task definitions for async processing
Future implementation for Redis RQ or Celery
"""
from typing import Dict, Optional


def predict_properties_task(smiles: str) -> Dict:
    """
    Background task for property prediction
    
    Args:
        smiles: SMILES string
        
    Returns:
        Prediction results
        
    Note:
        Placeholder for future async task implementation
    """
    # TODO: Implement with Redis RQ or Celery
    from backend.services.prediction_service import PredictionService
    return PredictionService.predict_properties(smiles)


def generate_molecule_task(prompt: str) -> Dict:
    """
    Background task for molecule generation
    
    Args:
        prompt: Generation prompt
        
    Returns:
        Generated molecule data
        
    Note:
        Placeholder for future async task implementation
    """
    # TODO: Implement with Redis RQ or Celery
    from backend.services.generation_service import GenerationService
    from backend.models.schemas.generation_schema import GenerateIn
    result = GenerationService.generate_molecule(GenerateIn(input=prompt))
    return result.model_dump()

