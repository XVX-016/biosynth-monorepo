"""
Energy calculation and geometry optimization API endpoints
"""
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List, Dict, Any, Optional
from chem.energy.forcefield import calculate_total_energy
from chem.energy.geometry import optimize_geometry
from chem.energy.dynamics import simulate_dynamics

router = APIRouter()

class MoleculeRequest(BaseModel):
    atoms: List[Dict[str, Any]]
    bonds: List[Dict[str, Any]]

class OptimizationRequest(MoleculeRequest):
    max_iterations: Optional[int] = 100
    convergence_threshold: Optional[float] = 0.001

class SimulationRequest(MoleculeRequest):
    temperature: Optional[float] = 300.0
    steps: Optional[int] = 100
    dt: Optional[float] = 0.001

@router.post("/calculate")
async def calculate_energy(req: MoleculeRequest):
    """Calculate molecular mechanics energy"""
    molecule = {
        "atoms": req.atoms,
        "bonds": req.bonds
    }
    try:
        result = calculate_total_energy(molecule)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/optimize")
async def optimize(req: OptimizationRequest):
    """Optimize molecular geometry"""
    molecule = {
        "atoms": req.atoms,
        "bonds": req.bonds
    }
    try:
        result = optimize_geometry(
            molecule,
            max_iterations=req.max_iterations or 100,
            convergence_threshold=req.convergence_threshold or 0.001
        )
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/simulate")
async def simulate(req: SimulationRequest):
    """Run molecular dynamics simulation"""
    molecule = {
        "atoms": req.atoms,
        "bonds": req.bonds
    }
    try:
        result = simulate_dynamics(
            molecule,
            temperature=req.temperature or 300.0,
            steps=req.steps or 100,
            dt=req.dt or 0.001
        )
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

