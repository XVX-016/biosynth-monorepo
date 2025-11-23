# backend/routes/library.py
from fastapi import APIRouter, Depends, HTTPException
from sqlmodel import Session
from typing import List
from backend.core.dependencies import get_db
from backend.services.molecule_service import MoleculeService
from backend.models.schemas.molecule_schema import MoleculeCreate, MoleculeResponse

router = APIRouter(prefix="/molecules", tags=["molecules"])

@router.post("/save", response_model=dict)
def save_molecule(payload: MoleculeCreate, db: Session = Depends(get_db)):
    """Save a molecule to the library"""
    molecule = MoleculeService.create_molecule(db, payload)
    return {
        "id": molecule.id,
        "name": molecule.name,
        "created_at": molecule.created_at.isoformat()
    }

@router.get("/list", response_model=List[dict])
def list_molecules(limit: int = 50, offset: int = 0, db: Session = Depends(get_db)):
    """List molecules with pagination"""
    molecules = MoleculeService.list_molecules(db, limit=limit, offset=offset)
    return [
        {
            "id": m.id,
            "name": m.name,
            "smiles": m.smiles,
            "properties": m.properties,
            "thumbnail_b64": m.thumbnail_b64,
            "molfile": m.molfile,
            "created_at": m.created_at.isoformat()
        }
        for m in molecules
    ]

@router.get("/{mol_id}", response_model=dict)
def get_molecule(mol_id: int, db: Session = Depends(get_db)):
    """Get a molecule by ID"""
    molecule = MoleculeService.get_molecule(db, mol_id)
    if not molecule:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return {
        "id": molecule.id,
        "name": molecule.name,
        "smiles": molecule.smiles,
        "json_graph": molecule.json_graph,
        "coords": molecule.coords,
        "properties": molecule.properties,
        "thumbnail_b64": molecule.thumbnail_b64,
        "molfile": molecule.molfile,
        "created_at": molecule.created_at.isoformat()
    }

@router.delete("/{mol_id}", response_model=dict)
def delete_molecule(mol_id: int, db: Session = Depends(get_db)):
    """Delete a molecule"""
    success = MoleculeService.delete_molecule(db, mol_id)
    if not success:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return {"status": "deleted", "id": mol_id}

