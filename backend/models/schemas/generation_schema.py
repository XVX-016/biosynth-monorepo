# backend/models/schemas/generation_schema.py
from pydantic import BaseModel
from typing import Optional


class GenerateIn(BaseModel):
    """Schema for molecule generation input"""
    input: str
    format: Optional[str] = "smiles"  # "smiles" or "text"


class GenerateOut(BaseModel):
    """Schema for molecule generation output"""
    smiles: Optional[str] = None
    json_graph: Optional[str] = None
    message: Optional[str] = None

