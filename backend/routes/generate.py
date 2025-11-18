"""
Molecule generation route handlers
"""
from fastapi import APIRouter, HTTPException
from backend.services.generation_service import GenerationService
from backend.models.schemas.generation_schema import GenerateIn, GenerateOut

router = APIRouter()


@router.post("/", response_model=GenerateOut)
def generate(payload: GenerateIn):
    """
    Generate molecule from prompt
    TODO: Implement transformer-based SMILES generation
    """
    result = GenerationService.generate_molecule(payload)
    return result

