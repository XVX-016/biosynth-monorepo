"""
Reaction simulation API endpoints
"""
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List, Dict, Any, Optional
from chem.reaction.engine import predict_reaction_products
from chem.reaction.mechanism import predict_mechanism
from chem.reaction.utils import calculate_reaction_energy

router = APIRouter()

class MoleculeRequest(BaseModel):
    atoms: List[Dict[str, Any]]
    bonds: List[Dict[str, Any]]

class ReactionRequest(BaseModel):
    reactants: List[MoleculeRequest]
    reaction_type: str
    conditions: Optional[Dict[str, Any]] = None

class MechanismRequest(BaseModel):
    reactants: List[MoleculeRequest]
    reaction_type: str
    max_steps: Optional[int] = 5

@router.post("/simulate")
async def simulate_reaction(req: ReactionRequest):
    """Simulate chemical reaction"""
    try:
        reactants = [
            {"atoms": r.atoms, "bonds": r.bonds}
            for r in req.reactants
        ]
        
        result = predict_reaction_products(
            reactants,
            req.reaction_type,
            req.conditions
        )
        
        # Calculate reaction energy if we have products
        if result["products"]:
            reaction_thermo = calculate_reaction_energy(reactants, result["products"])
            result["thermodynamics"] = reaction_thermo
        
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/mechanism")
async def predict_mechanism_endpoint(req: MechanismRequest):
    """Predict reaction mechanism"""
    try:
        reactants = [
            {"atoms": r.atoms, "bonds": r.bonds}
            for r in req.reactants
        ]
        
        result = predict_mechanism(
            reactants,
            req.reaction_type,
            req.max_steps or 5
        )
        
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

