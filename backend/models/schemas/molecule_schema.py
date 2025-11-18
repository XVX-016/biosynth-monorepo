# backend/models/schemas/molecule_schema.py
from pydantic import BaseModel
from typing import Optional, List
from datetime import datetime


class MoleculeCreate(BaseModel):
    """Schema for creating a molecule"""
    name: str
    smiles: Optional[str] = None
    json_graph: Optional[str] = None
    coords: Optional[str] = None
    properties: Optional[str] = None
    thumbnail_b64: Optional[str] = None


class MoleculeUpdate(BaseModel):
    """Schema for updating a molecule"""
    name: Optional[str] = None
    smiles: Optional[str] = None
    json_graph: Optional[str] = None
    coords: Optional[str] = None
    properties: Optional[str] = None
    thumbnail_b64: Optional[str] = None


class MoleculeResponse(BaseModel):
    """Schema for molecule response"""
    id: int
    name: str
    smiles: Optional[str] = None
    json_graph: Optional[str] = None
    coords: Optional[str] = None
    properties: Optional[str] = None
    thumbnail_b64: Optional[str] = None
    created_at: datetime
    owner_id: Optional[str] = None

    class Config:
        from_attributes = True


class ItemCreate(BaseModel):
    """Schema for creating an item"""
    name: str
    smiles: Optional[str] = None
    description: Optional[str] = None
    tags: Optional[List[str]] = None
    status: Optional[str] = 'in-stock'
    stock: Optional[int] = 0


class ItemUpdate(BaseModel):
    """Schema for updating an item"""
    name: Optional[str] = None
    smiles: Optional[str] = None
    description: Optional[str] = None
    tags: Optional[List[str]] = None
    status: Optional[str] = None
    stock: Optional[int] = None

