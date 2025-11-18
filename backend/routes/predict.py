"""
Prediction route handlers
"""
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List
from backend.services.prediction_service import PredictionService
from backend.models.schemas.prediction_schema import PredictIn, PredictOut
from backend.core.exceptions import InvalidSMILESError, PredictionError

router = APIRouter()


class BulkPredictIn(BaseModel):
    smiles_list: List[str]


class BulkPredictOut(BaseModel):
    results: List[dict]


@router.post("/", response_model=PredictOut)
def predict(payload: PredictIn):
    """
    Predict molecular properties from SMILES string
    """
    try:
        properties = PredictionService.predict_properties(payload.smiles, use_onnx=False)
        return PredictOut(properties=properties, smiles=payload.smiles)
    except InvalidSMILESError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except PredictionError as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/bulk", response_model=BulkPredictOut)
def predict_bulk(payload: BulkPredictIn):
    """
    Predict properties for multiple SMILES strings
    """
    results = []
    
    for smiles in payload.smiles_list:
        try:
            properties = PredictionService.predict_properties(smiles, use_onnx=False)
            results.append({
                "smiles": smiles,
                "error": None,
                "properties": properties
            })
        except (InvalidSMILESError, PredictionError) as e:
            results.append({
                "smiles": smiles,
                "error": str(e),
                "properties": None
            })
    
    return BulkPredictOut(results=results)

