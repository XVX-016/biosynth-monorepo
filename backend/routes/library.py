# backend/routes/library.py
from fastapi import APIRouter, Depends, HTTPException
from sqlmodel import select, Session
from typing import List
from backend.db import get_session
from backend.models_db import Molecule, MoleculeCreate

router = APIRouter(prefix="/molecules", tags=["molecules"])

@router.post("/save", response_model=dict)
def save_molecule(payload: MoleculeCreate, session: Session = Depends(get_session)):
    mol = Molecule(
        name=payload.name,
        smiles=payload.smiles,
        json_graph=payload.json_graph,
        coords=payload.coords,
        properties=payload.properties,
        thumbnail_b64=payload.thumbnail_b64,
    )
    session.add(mol)
    session.commit()
    session.refresh(mol)
    return {"id": mol.id, "name": mol.name, "created_at": mol.created_at.isoformat()}

@router.get("/list", response_model=List[dict])
def list_molecules(limit: int = 50, session: Session = Depends(get_session)):
    q = select(Molecule).order_by(Molecule.created_at.desc()).limit(limit)
    results = session.exec(q).all()
    out = []
    for m in results:
        out.append({
            "id": m.id,
            "name": m.name,
            "smiles": m.smiles,
            "properties": m.properties,
            "thumbnail_b64": m.thumbnail_b64,
            "created_at": m.created_at.isoformat()
        })
    return out

@router.get("/{mol_id}", response_model=dict)
def get_molecule(mol_id: int, session: Session = Depends(get_session)):
    m = session.get(Molecule, mol_id)
    if not m:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return {
        "id": m.id,
        "name": m.name,
        "smiles": m.smiles,
        "json_graph": m.json_graph,
        "coords": m.coords,
        "properties": m.properties,
        "thumbnail_b64": m.thumbnail_b64,
        "created_at": m.created_at.isoformat()
    }

@router.delete("/{mol_id}", response_model=dict)
def delete_molecule(mol_id: int, session: Session = Depends(get_session)):
    m = session.get(Molecule, mol_id)
    if not m:
        raise HTTPException(status_code=404, detail="Molecule not found")
    session.delete(m)
    session.commit()
    return {"status": "deleted", "id": mol_id}

