"""
Phase 8: Conformer Generation API routes

Provides endpoint for conformer generation.
"""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field
from typing import List, Dict
import logging
import sys
from pathlib import Path

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

from src.conformers import ConformerGenerator

logger = logging.getLogger(__name__)

router = APIRouter()

# Global instance
_conformer_generator = ConformerGenerator()


class ConformerGenerateRequest(BaseModel):
    """Request for conformer generation"""
    smiles: str = Field(..., description="SMILES string")
    n: int = Field(10, ge=1, le=100, description="Number of conformers")


@router.post("/generate")
async def generate_conformers(request: ConformerGenerateRequest):
    """
    Generate conformers for molecule.
    """
    try:
        conformers = _conformer_generator.generate_conformers(
            smiles=request.smiles,
            n=request.n
        )
        
        return {
            "conformers": conformers,
            "count": len(conformers),
            "smiles": request.smiles,
        }
    except Exception as e:
        logger.error(f"Conformer generation error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))

