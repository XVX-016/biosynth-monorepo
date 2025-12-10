from typing import Optional, Dict, Any, List
from pydantic import BaseModel
from datetime import datetime

class MoleculeBase(BaseModel):
    name: str
    smiles: Optional[str] = None
    json_graph: Optional[Dict[str, Any]] = None
    properties: Optional[Dict[str, Any]] = None
    thumbnail_b64: Optional[str] = None

class MoleculeCreate(MoleculeBase):
    pass

class MoleculeUpdate(BaseModel):
    name: Optional[str] = None
    smiles: Optional[str] = None
    json_graph: Optional[Dict[str, Any]] = None
    properties: Optional[Dict[str, Any]] = None
    thumbnail_b64: Optional[str] = None

class MoleculeResponse(MoleculeBase):
    id: int
    created_at: datetime
    updated_at: datetime
    
    class Config:
        from_attributes = True

# Aliases for legacy compatibility
ItemCreate = MoleculeCreate
ItemUpdate = MoleculeUpdate
