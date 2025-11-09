"""
Molecule generation route handlers
"""
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

router = APIRouter()


class GenerateIn(BaseModel):
    prompt: str


class GenerateOut(BaseModel):
    smiles: str


@router.post("/", response_model=GenerateOut)
def generate(payload: GenerateIn):
    """
    Generate molecule from prompt
    TODO: Implement transformer-based SMILES generation
    """
    # Placeholder: return simple carbon
    return GenerateOut(smiles="C")

