# backend/routes/library.py
from fastapi import APIRouter, Depends, HTTPException
from sqlmodel import Session
from typing import List
from pydantic import BaseModel
from backend.core.dependencies import get_db
from backend.services.molecule_service import MoleculeService
from backend.models.schemas.molecule_schema import MoleculeCreate, MoleculeResponse, MoleculeUpdate

router = APIRouter(prefix="/molecules", tags=["molecules"])

class MolfileUpdate(BaseModel):
    molfile: str

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
    try:
        molecules = MoleculeService.list_molecules(db, limit=limit, offset=offset)
        return [
            {
                "id": m.id,
                "name": m.name,
                "smiles": m.smiles,
                "formula": getattr(m, 'formula', None),  # Add formula if it exists
                "properties": m.properties,
                "thumbnail_b64": m.thumbnail_b64,
                "molfile": m.molfile,
                "created_at": m.created_at.isoformat()
            }
            for m in molecules
        ]
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to list molecules: {str(e)}")

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

@router.patch("/{mol_id}/molfile", response_model=dict)
def save_molfile(mol_id: int, payload: MolfileUpdate, db: Session = Depends(get_db)):
    """Save/update molfile for a molecule"""
    molecule = MoleculeService.get_molecule(db, mol_id)
    if not molecule:
        raise HTTPException(status_code=404, detail="Molecule not found")
    
    update_data = MoleculeUpdate(molfile=payload.molfile)
    updated = MoleculeService.update_molecule(db, mol_id, update_data)
    if not updated:
        raise HTTPException(status_code=500, detail="Failed to update molecule")
    
    return {"ok": True, "id": mol_id}

@router.post("/{mol_id}/export", response_model=dict)
def export_molecule(mol_id: int, db: Session = Depends(get_db)):
    """Export molecule for Lab (returns JSON graph)"""
    molecule = MoleculeService.get_molecule(db, mol_id)
    if not molecule:
        raise HTTPException(status_code=404, detail="Molecule not found")
    
    return {
        "id": molecule.id,
        "name": molecule.name,
        "formula": getattr(molecule, 'formula', None),
        "data": molecule.json_graph  # Lab expects this as 'data' or we just return graph
    }

@router.get("/{mol_id}/preview")
def preview_molecule(mol_id: int, db: Session = Depends(get_db)):
    """Generate or retrieve preview image"""
    molecule = MoleculeService.get_molecule(db, mol_id)
    if not molecule:
        raise HTTPException(status_code=404, detail="Molecule not found")

    # If we had a stored path, we'd use it. For now, on-the-fly render
    if molecule.json_graph:
        from backend.utils import preview as utils_preview
        from fastapi.responses import StreamingResponse
        from io import BytesIO
        
        png_bytes = utils_preview.render_preview_png(molecule.json_graph)
        return StreamingResponse(BytesIO(png_bytes), media_type="image/png")
    
    raise HTTPException(status_code=404, detail="No structure data for preview")
