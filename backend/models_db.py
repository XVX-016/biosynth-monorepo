# backend/models_db.py
from sqlmodel import SQLModel, Field
from typing import Optional, List
from datetime import datetime
from pydantic import BaseModel
import json

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

class Item(SQLModel, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    name: str
    smiles: Optional[str] = None
    description: Optional[str] = None
    tags: Optional[str] = None  # JSON string of tags array
    status: str = Field(default='in-stock')  # 'in-stock', 'sold-out', 'archived'
    stock: Optional[int] = Field(default=0)
    structure_file_url: Optional[str] = None
    created_at: datetime = Field(default_factory=datetime.utcnow)
    updated_at: datetime = Field(default_factory=datetime.utcnow)
    
    def get_tags(self) -> List[str]:
        """Get tags as list"""
        if not self.tags:
            return []
        try:
            return json.loads(self.tags)
        except:
            return []
    
    def set_tags(self, tags: List[str]):
        """Set tags from list"""
        self.tags = json.dumps(tags) if tags else None

class ItemCreate(BaseModel):
    name: str
    smiles: Optional[str] = None
    description: Optional[str] = None
    tags: Optional[List[str]] = None
    status: Optional[str] = 'in-stock'
    stock: Optional[int] = 0

class ItemUpdate(BaseModel):
    name: Optional[str] = None
    smiles: Optional[str] = None
    description: Optional[str] = None
    tags: Optional[List[str]] = None
    status: Optional[str] = None
    stock: Optional[int] = None

