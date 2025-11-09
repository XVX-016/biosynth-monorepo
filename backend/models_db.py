# backend/models_db.py
from sqlmodel import SQLModel, Field
from typing import Optional
from datetime import datetime
from pydantic import BaseModel

class Molecule(SQLModel, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    name: str
    smiles: Optional[str] = None
    json_graph: Optional[str] = None   # store JSON stringified engine graph
    coords: Optional[str] = None       # optional stored coordinates (JSON)
    properties: Optional[str] = None   # JSON string of properties
    thumbnail_b64: Optional[str] = None
    created_at: datetime = Field(default_factory=datetime.utcnow)
    owner_id: Optional[str] = None     # optional user id if auth added later

class MoleculeCreate(BaseModel):
    name: str
    smiles: Optional[str] = None
    json_graph: Optional[str] = None
    coords: Optional[str] = None
    properties: Optional[str] = None
    thumbnail_b64: Optional[str] = None

