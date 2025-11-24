"""
Spectroscopy API endpoints
"""
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List, Dict, Any
from chem.spectroscopy.ir import predict_ir_spectrum
from chem.spectroscopy.nmr import predict_nmr_spectrum
from chem.spectroscopy.mass import predict_mass_spectrum

router = APIRouter()

class MoleculeRequest(BaseModel):
    atoms: List[Dict[str, Any]]
    bonds: List[Dict[str, Any]]

@router.post("/ir")
async def predict_ir(req: MoleculeRequest):
    """Predict IR spectrum"""
    molecule = {
        "atoms": req.atoms,
        "bonds": req.bonds
    }
    return predict_ir_spectrum(molecule)

@router.post("/nmr")
async def predict_nmr(req: MoleculeRequest):
    """Predict NMR spectrum (1H and 13C)"""
    molecule = {
        "atoms": req.atoms,
        "bonds": req.bonds
    }
    return predict_nmr_spectrum(molecule)

@router.post("/mass")
async def predict_mass(req: MoleculeRequest):
    """Predict Mass Spectrometry spectrum"""
    molecule = {
        "atoms": req.atoms,
        "bonds": req.bonds
    }
    return predict_mass_spectrum(molecule)

