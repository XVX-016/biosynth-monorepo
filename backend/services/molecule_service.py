# backend/services/molecule_service.py
"""
Service layer for molecule CRUD operations
"""
from typing import List, Optional
from sqlmodel import Session, select
from backend.models.db.molecule import Molecule
from backend.models.schemas.molecule_schema import MoleculeCreate, MoleculeUpdate


class MoleculeService:
    """Service for molecule management"""
    
    @staticmethod
    def create_molecule(db: Session, molecule_data: MoleculeCreate) -> Molecule:
        """
        Create a new molecule
        
        Args:
            db: Database session
            molecule_data: Molecule creation data
            
        Returns:
            Created molecule
        """
        molecule = Molecule(**molecule_data.model_dump())
        db.add(molecule)
        db.commit()
        db.refresh(molecule)
        return molecule
    
    @staticmethod
    def get_molecule(db: Session, molecule_id: int) -> Optional[Molecule]:
        """
        Get molecule by ID
        
        Args:
            db: Database session
            molecule_id: Molecule ID
            
        Returns:
            Molecule or None if not found
        """
        return db.get(Molecule, molecule_id)
    
    @staticmethod
    def list_molecules(db: Session, limit: int = 100, offset: int = 0) -> List[Molecule]:
        """
        List molecules with pagination
        
        Args:
            db: Database session
            limit: Maximum number of molecules to return
            offset: Number of molecules to skip
            
        Returns:
            List of molecules
        """
        statement = select(Molecule).offset(offset).limit(limit).order_by(Molecule.created_at.desc())
        return list(db.exec(statement))
    
    @staticmethod
    def update_molecule(
        db: Session,
        molecule_id: int,
        molecule_data: MoleculeUpdate
    ) -> Optional[Molecule]:
        """
        Update a molecule
        
        Args:
            db: Database session
            molecule_id: Molecule ID
            molecule_data: Update data
            
        Returns:
            Updated molecule or None if not found
        """
        molecule = db.get(Molecule, molecule_id)
        if molecule is None:
            return None
        
        update_data = molecule_data.model_dump(exclude_unset=True)
        for key, value in update_data.items():
            setattr(molecule, key, value)
        
        db.add(molecule)
        db.commit()
        db.refresh(molecule)
        return molecule
    
    @staticmethod
    def delete_molecule(db: Session, molecule_id: int) -> bool:
        """
        Delete a molecule
        
        Args:
            db: Database session
            molecule_id: Molecule ID
            
        Returns:
            True if deleted, False if not found
        """
        molecule = db.get(Molecule, molecule_id)
        if molecule is None:
            return False
        
        db.delete(molecule)
        db.commit()
        return True

