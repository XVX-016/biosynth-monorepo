"""
Phase 8: Quantum Chemistry and Molecular Dynamics API routes

Provides endpoints for QM calculations, MD simulations, and conformer generation.
"""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field
from typing import List, Optional, Dict, Any
import logging
import sys
from pathlib import Path

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

from src.qm import QMEngine, Psi4Wrapper, XTBWrapper
from src.md import MDEngine
from src.conformers import ConformerGenerator, ETKDGGenerator

logger = logging.getLogger(__name__)

router = APIRouter()

# Global instances
_qm_engine = QMEngine()
_md_engine = MDEngine()
_conformer_generator = ConformerGenerator()


# Request models
class QMEnergyRequest(BaseModel):
    smiles: str
    coordinates: Optional[List[List[float]]] = None
    method: str = Field("B3LYP", description="QM method")
    basis: str = Field("6-31G", description="Basis set")


class QMOptimizeRequest(BaseModel):
    smiles: str
    initial_coordinates: Optional[List[List[float]]] = None
    method: str = Field("B3LYP", description="QM method")
    basis: str = Field("6-31G", description="Basis set")


class MDSimulateRequest(BaseModel):
    smiles: str
    coordinates: Optional[List[List[float]]] = None
    steps: int = Field(1000, ge=1, le=100000, description="Number of MD steps")
    timestep: float = Field(0.001, description="Time step in picoseconds")
    temperature: float = Field(300.0, description="Temperature in Kelvin")
    forcefield: str = Field("UFF", description="Force field name")


class ConformerGenerateRequest(BaseModel):
    smiles: str
    n: int = Field(10, ge=1, le=100, description="Number of conformers")
    initial_coordinates: Optional[List[List[float]]] = None


@router.post("/qm/energy")
async def compute_qm_energy(request: QMEnergyRequest):
    """
    Compute quantum mechanical energy.
    """
    try:
        result = await _qm_engine.compute_energy(
            smiles=request.smiles,
            coordinates=request.coordinates,
            method=request.method,
            basis=request.basis
        )
        
        return {
            "energy": result.energy,
            "coordinates": result.coordinates,
            "metadata": result.metadata,
        }
    except Exception as e:
        logger.error(f"QM energy error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/qm/optimize")
async def optimize_geometry(request: QMOptimizeRequest):
    """
    Optimize molecular geometry using QM.
    """
    try:
        result = await _qm_engine.optimize_geometry(
            smiles=request.smiles,
            initial_coordinates=request.initial_coordinates,
            method=request.method,
            basis=request.basis
        )
        
        return {
            "energy": result.energy,
            "coordinates": result.coordinates,
            "forces": result.forces,
            "metadata": result.metadata,
        }
    except Exception as e:
        logger.error(f"QM optimize error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/md/simulate")
async def simulate_md(request: MDSimulateRequest):
    """
    Run molecular dynamics simulation.
    """
    try:
        result = await _md_engine.simulate(
            smiles=request.smiles,
            coordinates=request.coordinates,
            steps=request.steps,
            timestep=request.timestep,
            temperature=request.temperature,
            forcefield=request.forcefield
        )
        
        # Convert trajectory frames to dicts
        trajectory = [
            {
                "step": frame.step,
                "time": frame.time,
                "coordinates": frame.coordinates,
                "velocities": frame.velocities,
                "energy": frame.energy,
                "temperature": frame.temperature,
            }
            for frame in result.trajectory
        ]
        
        return {
            "trajectory": trajectory,
            "final_coordinates": result.final_coordinates,
            "total_energy": result.total_energy,
            "metadata": result.metadata,
        }
    except Exception as e:
        logger.error(f"MD simulate error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/conformers/generate")
async def generate_conformers(request: ConformerGenerateRequest):
    """
    Generate conformers for molecule.
    """
    try:
        conformers = _conformer_generator.generate_conformers(
            smiles=request.smiles,
            n=request.n,
            initial_coordinates=request.initial_coordinates
        )
        
        return {
            "conformers": conformers,
            "count": len(conformers),
            "smiles": request.smiles,
        }
    except Exception as e:
        logger.error(f"Conformer generation error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))

