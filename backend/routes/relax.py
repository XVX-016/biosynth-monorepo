# backend/routes/relax.py
from fastapi import APIRouter
from pydantic import BaseModel
from typing import List

router = APIRouter()


class AtomIn(BaseModel):
    id: str
    element: str
    position: List[float]


class BondIn(BaseModel):
    id: str
    atom1: str
    atom2: str
    order: int


class MolIn(BaseModel):
    atoms: List[AtomIn]
    bonds: List[BondIn]


@router.post("/relax")
def relax(mol: MolIn):
    """
    Relax molecule geometry.
    
    Placeholder: immediate echo. Replace with RDKit ETKDG + MMFF optimization later.
    For now, returns atoms with x/y/z coordinates (same as input).
    
    TODO: Implement RDKit:
    1. Convert incoming molecule to RDKit Mol object
    2. Use ETKDG to generate 3D coordinates
    3. Optimize with MMFF
    4. Return optimized coordinates + energy
    """
    # Placeholder: return same coordinates
    out_atoms = []
    for a in mol.atoms:
        out_atoms.append({
            "id": a.id,
            "x": a.position[0],
            "y": a.position[1],
            "z": a.position[2]
        })
    return {"atoms": out_atoms, "energy": 0.0}

