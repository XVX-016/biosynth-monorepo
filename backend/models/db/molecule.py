# backend/models/db/molecule.py
from sqlmodel import SQLModel, Field
from typing import Optional, List
from datetime import datetime
import json


class Molecule(SQLModel, table=True):
    """Database model for molecules"""
    id: Optional[int] = Field(default=None, primary_key=True)
    name: str
    smiles: Optional[str] = None
    json_graph: Optional[str] = None   # store JSON stringified engine graph
    coords: Optional[str] = None       # optional stored coordinates (JSON)
    properties: Optional[str] = None   # JSON string of properties
    thumbnail_b64: Optional[str] = None
    created_at: datetime = Field(default_factory=datetime.utcnow)
    owner_id: Optional[str] = None     # optional user id if auth added later


class Item(SQLModel, table=True):
    """Database model for items (admin/inventory)"""
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

