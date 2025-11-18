"""
Admin routes for item management
"""
import os
import json
from fastapi import APIRouter, Depends, HTTPException, UploadFile, File, Form
from sqlmodel import select, Session
from typing import List, Optional
from backend.core.dependencies import get_db
from backend.models.db.molecule import Item
from backend.models.schemas.molecule_schema import ItemCreate, ItemUpdate

router = APIRouter(prefix="/api/v1/admin/items", tags=["admin"])

# File upload directory
UPLOAD_DIR = os.path.join(os.path.dirname(__file__), "../../uploads")
os.makedirs(UPLOAD_DIR, exist_ok=True)


def item_to_dict(item: Item) -> dict:
    """Convert Item model to dictionary"""
    return {
        "id": item.id,
        "name": item.name,
        "smiles": item.smiles,
        "description": item.description,
        "tags": item.get_tags(),
        "status": item.status,
        "stock": item.stock,
        "structure_file_url": item.structure_file_url,
        "created_at": item.created_at.isoformat(),
        "updated_at": item.updated_at.isoformat(),
    }


@router.get("", response_model=List[dict])
def list_items(session: Session = Depends(get_db)):
    """List all items"""
    query = select(Item).order_by(Item.created_at.desc())
    items = session.exec(query).all()
    return [item_to_dict(item) for item in items]


@router.get("/{item_id}", response_model=dict)
def get_item(item_id: int, session: Session = Depends(get_db)):
    """Get item by ID"""
    item = session.get(Item, item_id)
    if not item:
        raise HTTPException(status_code=404, detail="Item not found")
    return item_to_dict(item)


@router.post("", response_model=dict)
def create_item(
    name: str = Form(...),
    smiles: Optional[str] = Form(None),
    description: Optional[str] = Form(None),
    tags: Optional[str] = Form(None),
    status: Optional[str] = Form("in-stock"),
    stock: Optional[int] = Form(0),
    structure_file: Optional[UploadFile] = File(None),
    session: Session = Depends(get_db),
):
    """Create new item"""
    # Parse tags
    tags_list = []
    if tags:
        try:
            tags_list = json.loads(tags)
        except:
            pass
    
    # Handle file upload
    structure_file_url = None
    if structure_file:
        # Save file
        file_path = os.path.join(UPLOAD_DIR, f"{name}_{structure_file.filename}")
        with open(file_path, "wb") as f:
            content = structure_file.file.read()
            f.write(content)
        structure_file_url = f"/uploads/{os.path.basename(file_path)}"
    
    # Create item
    item = Item(
        name=name,
        smiles=smiles,
        description=description,
        status=status or "in-stock",
        stock=stock or 0,
        structure_file_url=structure_file_url,
    )
    item.set_tags(tags_list)
    
    session.add(item)
    session.commit()
    session.refresh(item)
    
    return item_to_dict(item)


@router.put("/{item_id}", response_model=dict)
def update_item(
    item_id: int,
    name: Optional[str] = Form(None),
    smiles: Optional[str] = Form(None),
    description: Optional[str] = Form(None),
    tags: Optional[str] = Form(None),
    status: Optional[str] = Form(None),
    stock: Optional[int] = Form(None),
    structure_file: Optional[UploadFile] = File(None),
    session: Session = Depends(get_db),
):
    """Update item"""
    item = session.get(Item, item_id)
    if not item:
        raise HTTPException(status_code=404, detail="Item not found")
    
    # Update fields
    if name is not None:
        item.name = name
    if smiles is not None:
        item.smiles = smiles
    if description is not None:
        item.description = description
    if status is not None:
        item.status = status
    if stock is not None:
        item.stock = stock
    
    # Update tags
    if tags is not None:
        try:
            tags_list = json.loads(tags)
            item.set_tags(tags_list)
        except:
            pass
    
    # Handle file upload
    if structure_file:
        # Delete old file if exists
        if item.structure_file_url:
            old_path = os.path.join(UPLOAD_DIR, os.path.basename(item.structure_file_url))
            if os.path.exists(old_path):
                os.remove(old_path)
        
        # Save new file
        file_path = os.path.join(UPLOAD_DIR, f"{item.name}_{structure_file.filename}")
        with open(file_path, "wb") as f:
            content = structure_file.file.read()
            f.write(content)
        item.structure_file_url = f"/uploads/{os.path.basename(file_path)}"
    
    from datetime import datetime
    item.updated_at = datetime.utcnow()
    
    session.add(item)
    session.commit()
    session.refresh(item)
    
    return item_to_dict(item)


@router.delete("/{item_id}", response_model=dict)
def delete_item(item_id: int, session: Session = Depends(get_db)):
    """Delete item"""
    item = session.get(Item, item_id)
    if not item:
        raise HTTPException(status_code=404, detail="Item not found")
    
    # Delete associated file
    if item.structure_file_url:
        file_path = os.path.join(UPLOAD_DIR, os.path.basename(item.structure_file_url))
        if os.path.exists(file_path):
            os.remove(file_path)
    
    session.delete(item)
    session.commit()
    
    return {"status": "deleted", "id": item_id}

