"""
Retrosynthesis planning API endpoints
"""
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List, Dict, Any, Optional
from chem.retrosynthesis.planner import plan_retrosynthesis
from chem.retrosynthesis.scoring import rank_pathways

router = APIRouter()

class MoleculeRequest(BaseModel):
    atoms: List[Dict[str, Any]]
    bonds: List[Dict[str, Any]]

class RetrosynthesisRequest(MoleculeRequest):
    max_steps: Optional[int] = 5
    max_pathways: Optional[int] = 10

@router.post("/plan")
async def plan_retrosynthesis_route(req: RetrosynthesisRequest):
    """Plan retrosynthesis pathways for target molecule"""
    target_molecule = {
        "atoms": req.atoms,
        "bonds": req.bonds
    }
    
    try:
        # Generate pathways
        pathways = plan_retrosynthesis(
            target_molecule,
            max_steps=req.max_steps or 5,
            max_pathways=req.max_pathways or 10
        )
        
        # Score and rank pathways
        ranked_pathways = rank_pathways(pathways)
        
        return {
            "target": target_molecule,
            "pathways": ranked_pathways,
            "count": len(ranked_pathways)
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

