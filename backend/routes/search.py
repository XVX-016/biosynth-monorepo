"""
Search API endpoints for substructure and SMARTS matching
"""
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List, Dict, Any, Optional
from chem.search.substructure import find_substructure
from chem.search.smarts import match_smarts
from chem.search.indexer import get_index

router = APIRouter()

class MoleculeQuery(BaseModel):
    atoms: List[Dict[str, Any]]
    bonds: List[Dict[str, Any]]

class SMARTSQuery(BaseModel):
    pattern: str

class SubstructureSearchRequest(BaseModel):
    query: MoleculeQuery

class SMARTSSearchRequest(BaseModel):
    pattern: str

@router.post("/substructure")
async def search_substructure(req: SubstructureSearchRequest):
    """
    Search for molecules containing the query substructure
    """
    query_mol = {
        "atoms": req.query.atoms,
        "bonds": req.query.bonds
    }
    
    # Get index
    index = get_index()
    candidates = index.search_candidates(query_mol)
    
    # Match against candidates
    matches = []
    for mol_id in candidates:
        target_mol = index.get_molecule(mol_id)
        if not target_mol:
            continue
        
        is_match, mapping = find_substructure(query_mol, {
            "atoms": target_mol["atoms"],
            "bonds": target_mol["bonds"]
        })
        
        if is_match:
            matches.append({
                "molecule_id": mol_id,
                "matched_atoms": [{"query": p, "target": t} for p, t in (mapping or [])]
            })
    
    return {"matches": matches, "count": len(matches)}

@router.post("/smarts")
async def search_smarts(req: SMARTSSearchRequest):
    """
    Search using SMARTS pattern
    """
    pattern = req.pattern
    
    # Get index
    index = get_index()
    all_molecules = list(index.molecules.keys())
    
    matches = []
    for mol_id in all_molecules:
        target_mol = index.get_molecule(mol_id)
        if not target_mol:
            continue
        
        is_match, matched_atoms = match_smarts(pattern, {
            "atoms": target_mol["atoms"],
            "bonds": target_mol["bonds"]
        })
        
        if is_match:
            matches.append({
                "molecule_id": mol_id,
                "matched_atoms": matched_atoms or []
            })
    
    return {"matches": matches, "count": len(matches)}

