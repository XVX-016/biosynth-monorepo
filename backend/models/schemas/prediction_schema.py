# backend/models/schemas/prediction_schema.py
from pydantic import BaseModel
from typing import Optional


class PredictIn(BaseModel):
    """Schema for prediction input"""
    smiles: str


class PredictOut(BaseModel):
    """Schema for prediction output"""
    properties: dict
    smiles: Optional[str] = None

