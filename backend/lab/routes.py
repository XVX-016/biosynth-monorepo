from fastapi import APIRouter
from pydantic import BaseModel
from typing import List, Tuple, Dict, Optional

router = APIRouter()

class AtomIn(BaseModel):
    id: str
    element: str
    position: Tuple[float,float,float]

class PredictRequest(BaseModel):
    atoms: List[AtomIn]

class BondOut(BaseModel):
    id: str
    a: str
    b: str
    order: int

class PredictResponse(BaseModel):
    bonds: List[BondOut]

class ExportRequest(BaseModel):
    atoms: List[Dict]
    bonds: List[Dict]
    format: str  # 'molforge', 'mol', or 'pdb'

class ExportResponse(BaseModel):
    data: str
    format: str

# Import bond engine and serializer
from chem.engine.bond_engine import BondEngine
from chem.core.serializer import MoleculeSerializer
from chem.core.models import Atom, Bond

bond_engine = BondEngine()
serializer = MoleculeSerializer()

@router.post('/predict-bonds', response_model=PredictResponse)
async def predict_bonds(req: PredictRequest):
    """Predict bonds using the full bond engine pipeline."""
    # Convert request to internal format
    atoms = {}
    positions = {}
    
    for a in req.atoms:
        atom_id = int(a.id) if a.id.isdigit() else hash(a.id) % 100000
        atoms[atom_id] = Atom(id=atom_id, element=a.element)
        positions[atom_id] = a.position
    
    # Run full bond prediction pipeline
    bonds = bond_engine.predict_all_bonds(atoms, positions)
    
    # Convert to response format
    bond_out = [
        BondOut(
            id=f"b_{b.atom_a}_{b.atom_b}",
            a=str(b.atom_a),
            b=str(b.atom_b),
            order=b.order
        )
        for b in bonds
    ]
    
    return PredictResponse(bonds=bond_out)

@router.post('/export-molecule', response_model=ExportResponse)
async def export_molecule(req: ExportRequest):
    """Export molecule to various formats."""
    # Convert to internal format
    atoms = {}
    positions = {}
    
    for a in req.atoms:
        atom_id = int(a['id']) if str(a['id']).isdigit() else hash(str(a['id'])) % 100000
        atoms[atom_id] = Atom(
            id=atom_id,
            element=a['element'],
            charge=a.get('charge', 0)
        )
        positions[atom_id] = (
            a.get('x', 0.0),
            a.get('y', 0.0),
            a.get('z', 0.0)
        )
    
    bonds = [
        Bond(
            atom_a=int(b['a']) if str(b['a']).isdigit() else hash(str(b['a'])) % 100000,
            atom_b=int(b['b']) if str(b['b']).isdigit() else hash(str(b['b'])) % 100000,
            order=b.get('order', 1)
        )
        for b in req.bonds
    ]
    
    # Export based on format
    if req.format == 'molforge':
        data = serializer.to_molforge(atoms, bonds, positions)
    elif req.format == 'mol':
        data = serializer.to_mol(atoms, bonds, positions)
    elif req.format == 'pdb':
        data = serializer.to_pdb(atoms, positions)
    else:
        data = serializer.to_molforge(atoms, bonds, positions)
    
    return ExportResponse(data=data, format=req.format)

@router.post('/parse-molecule')
async def parse_molecule():
    # TODO: Implement parse-molecule
    return {"status": "not implemented"}

@router.post('/save-session')
async def save_session():
    # TODO: Implement save-session
    return {"status": "not implemented"}

@router.post('/run-ml-prediction')
async def run_ml_prediction():
    # TODO: Implement run-ml-prediction
    return {"status": "not implemented"}

