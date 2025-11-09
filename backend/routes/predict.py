"""
Prediction route handlers
"""
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from backend.models.predictor import get_predictor
from backend.utils.featurizer import featurize_smiles, validate_smiles

router = APIRouter()


class PredictIn(BaseModel):
    smiles: str


class PredictOut(BaseModel):
    properties: dict


@router.post("/", response_model=PredictOut)
def predict(payload: PredictIn):
    """
    Predict molecular properties from SMILES string
    """
    if not validate_smiles(payload.smiles):
        raise HTTPException(status_code=400, detail="Invalid SMILES string")
    
    features = featurize_smiles(payload.smiles)
    if features is None:
        raise HTTPException(status_code=400, detail="Failed to featurize SMILES")
    
    predictor = get_predictor()
    properties = predictor.predict(features)
    
    return PredictOut(properties=properties)

