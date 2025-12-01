"""
QM Engine Interfaces - Protocols for quantum chemistry engines

Defines interfaces that real QM engines (Psi4, XTB, etc.) should implement.
"""

from typing import Protocol, Dict, List, Optional, Any
from dataclasses import dataclass
import numpy as np


@dataclass
class QMResult:
    """Result from QM calculation"""
    energy: float
    coordinates: Optional[List[List[float]]] = None  # [[x, y, z], ...] per atom
    forces: Optional[List[List[float]]] = None
    dipole: Optional[List[float]] = None
    charges: Optional[List[float]] = None
    metadata: Dict[str, Any] = None
    
    def __post_init__(self):
        if self.metadata is None:
            self.metadata = {}


class QMEngineProtocol(Protocol):
    """
    Protocol for quantum chemistry engines.
    
    Real engines (Psi4, XTB, etc.) should implement this interface.
    """
    
    async def compute_energy(
        self,
        smiles: str,
        coordinates: Optional[List[List[float]]] = None,
        method: str = "B3LYP",
        basis: str = "6-31G"
    ) -> QMResult:
        """
        Compute electronic energy.
        
        Args:
            smiles: SMILES string
            coordinates: Optional 3D coordinates
            method: QM method (e.g., "B3LYP", "HF", "MP2")
            basis: Basis set (e.g., "6-31G", "cc-pVDZ")
        
        Returns:
            QMResult with energy
        """
        ...
    
    async def optimize_geometry(
        self,
        smiles: str,
        initial_coordinates: Optional[List[List[float]]] = None,
        method: str = "B3LYP",
        basis: str = "6-31G"
    ) -> QMResult:
        """
        Optimize molecular geometry.
        
        Args:
            smiles: SMILES string
            initial_coordinates: Optional starting coordinates
            method: QM method
            basis: Basis set
        
        Returns:
            QMResult with optimized coordinates and energy
        """
        ...
    
    async def compute_properties(
        self,
        smiles: str,
        coordinates: Optional[List[List[float]]] = None,
        method: str = "B3LYP",
        basis: str = "6-31G"
    ) -> QMResult:
        """
        Compute molecular properties (dipole, charges, etc.).
        
        Args:
            smiles: SMILES string
            coordinates: Optional 3D coordinates
            method: QM method
            basis: Basis set
        
        Returns:
            QMResult with properties
        """
        ...

