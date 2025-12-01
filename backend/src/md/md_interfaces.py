"""
MD Engine Interfaces - Protocols for molecular dynamics engines

Defines interfaces that real MD engines (OpenMM, GROMACS, etc.) should implement.
"""

from typing import Protocol, List, Optional, Dict, Any
from dataclasses import dataclass
import numpy as np


@dataclass
class TrajectoryFrame:
    """Single frame from MD trajectory"""
    step: int
    time: float  # picoseconds
    coordinates: List[List[float]]  # [[x, y, z], ...] per atom
    velocities: Optional[List[List[float]]] = None
    energy: Optional[float] = None
    temperature: Optional[float] = None
    pressure: Optional[float] = None


@dataclass
class MDResult:
    """Result from MD simulation"""
    trajectory: List[TrajectoryFrame]
    final_coordinates: List[List[float]]
    total_energy: float
    metadata: Dict[str, Any] = None
    
    def __post_init__(self):
        if self.metadata is None:
            self.metadata = {}


class MDEngineProtocol(Protocol):
    """
    Protocol for molecular dynamics engines.
    
    Real engines (OpenMM, GROMACS, etc.) should implement this interface.
    """
    
    async def simulate(
        self,
        smiles: str,
        coordinates: Optional[List[List[float]]] = None,
        steps: int = 1000,
        timestep: float = 0.001,  # picoseconds
        temperature: float = 300.0,  # Kelvin
        forcefield: str = "UFF"
    ) -> MDResult:
        """
        Run MD simulation.
        
        Args:
            smiles: SMILES string
            coordinates: Optional initial 3D coordinates
            steps: Number of simulation steps
            timestep: Time step in picoseconds
            temperature: Target temperature in Kelvin
            forcefield: Force field name
        
        Returns:
            MDResult with trajectory
        """
        ...

