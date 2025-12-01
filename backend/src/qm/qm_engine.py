"""
QM Engine - Mock quantum chemistry engine

Provides mock implementations of QM calculations.
All results are fake but follow realistic patterns.
"""

from typing import List, Optional, Dict, Any
import numpy as np
import random
import logging

from .qm_interfaces import QMEngineProtocol, QMResult

logger = logging.getLogger(__name__)


class QMEngine:
    """
    Mock quantum chemistry engine.
    
    Returns fake but realistic-looking results.
    """
    
    def __init__(self):
        self.random_seed = None
    
    def _set_seed(self, smiles: str):
        """Set random seed based on SMILES for reproducibility."""
        self.random_seed = hash(smiles) % 1000000
        random.seed(self.random_seed)
        np.random.seed(self.random_seed)
    
    def _estimate_atom_count(self, smiles: str) -> int:
        """Estimate number of atoms from SMILES."""
        # Simple heuristic: count capital letters (rough atom count)
        return sum(1 for c in smiles if c.isupper() and c.isalpha())
    
    def _generate_fake_coordinates(
        self,
        n_atoms: int,
        smiles: str
    ) -> List[List[float]]:
        """
        Generate fake 3D coordinates.
        
        Creates a roughly spherical arrangement.
        """
        coords = []
        radius = 2.0  # Angstrom
        
        for i in range(n_atoms):
            # Generate points on sphere
            theta = 2 * np.pi * i / n_atoms if n_atoms > 0 else 0
            phi = np.pi * random.random()
            
            x = radius * np.sin(phi) * np.cos(theta) + random.gauss(0, 0.5)
            y = radius * np.sin(phi) * np.sin(theta) + random.gauss(0, 0.5)
            z = radius * np.cos(phi) + random.gauss(0, 0.5)
            
            coords.append([float(x), float(y), float(z)])
        
        return coords
    
    def _compute_fake_energy(
        self,
        smiles: str,
        method: str = "B3LYP"
    ) -> float:
        """
        Compute fake energy value.
        
        Uses heuristics to generate realistic-looking energies.
        """
        n_atoms = self._estimate_atom_count(smiles)
        
        # Base energy per atom (Hartree)
        base_energy_per_atom = -37.0  # Typical for organic molecules
        
        # Method corrections
        method_corrections = {
            "HF": 0.0,
            "B3LYP": -0.1,  # DFT typically lower
            "MP2": -0.15,
            "CCSD": -0.2,
        }
        correction = method_corrections.get(method, 0.0)
        
        # Add some randomness
        noise = random.gauss(0, 0.5)
        
        total_energy = n_atoms * base_energy_per_atom + correction * n_atoms + noise
        
        return round(total_energy, 6)
    
    async def compute_energy(
        self,
        smiles: str,
        coordinates: Optional[List[List[float]]] = None,
        method: str = "B3LYP",
        basis: str = "6-31G"
    ) -> QMResult:
        """
        Compute electronic energy (mock).
        
        Args:
            smiles: SMILES string
            coordinates: Optional 3D coordinates
            method: QM method
            basis: Basis set
        
        Returns:
            QMResult with energy
        """
        self._set_seed(smiles)
        
        energy = self._compute_fake_energy(smiles, method)
        
        # Generate coordinates if not provided
        if coordinates is None:
            n_atoms = self._estimate_atom_count(smiles)
            coordinates = self._generate_fake_coordinates(n_atoms, smiles)
        
        return QMResult(
            energy=energy,
            coordinates=coordinates,
            metadata={
                'method': method,
                'basis': basis,
                'smiles': smiles,
                'mock': True,
            }
        )
    
    async def optimize_geometry(
        self,
        smiles: str,
        initial_coordinates: Optional[List[List[float]]] = None,
        method: str = "B3LYP",
        basis: str = "6-31G"
    ) -> QMResult:
        """
        Optimize molecular geometry (mock).
        
        Args:
            smiles: SMILES string
            initial_coordinates: Optional starting coordinates
            method: QM method
            basis: Basis set
        
        Returns:
            QMResult with optimized coordinates and energy
        """
        self._set_seed(smiles)
        
        # Generate optimized coordinates (slightly more compact)
        n_atoms = self._estimate_atom_count(smiles)
        
        if initial_coordinates is None:
            coords = self._generate_fake_coordinates(n_atoms, smiles)
        else:
            # "Optimize" by slightly adjusting coordinates
            coords = []
            for x, y, z in initial_coordinates:
                # Move atoms slightly closer together (optimization effect)
                scale = 0.95
                coords.append([
                    float(x * scale + random.gauss(0, 0.1)),
                    float(y * scale + random.gauss(0, 0.1)),
                    float(z * scale + random.gauss(0, 0.1)),
                ])
        
        # Compute energy (should be lower after optimization)
        energy = self._compute_fake_energy(smiles, method)
        energy = energy - abs(random.gauss(0.1, 0.05))  # Lower energy after optimization
        
        # Generate fake forces (should be small after optimization)
        forces = [
            [random.gauss(0, 0.01), random.gauss(0, 0.01), random.gauss(0, 0.01)]
            for _ in range(n_atoms)
        ]
        
        return QMResult(
            energy=round(energy, 6),
            coordinates=coords,
            forces=forces,
            metadata={
                'method': method,
                'basis': basis,
                'smiles': smiles,
                'optimized': True,
                'mock': True,
            }
        )
    
    async def compute_properties(
        self,
        smiles: str,
        coordinates: Optional[List[List[float]]] = None,
        method: str = "B3LYP",
        basis: str = "6-31G"
    ) -> QMResult:
        """
        Compute molecular properties (mock).
        
        Args:
            smiles: SMILES string
            coordinates: Optional 3D coordinates
            method: QM method
            basis: Basis set
        
        Returns:
            QMResult with properties
        """
        self._set_seed(smiles)
        
        n_atoms = self._estimate_atom_count(smiles)
        
        if coordinates is None:
            coordinates = self._generate_fake_coordinates(n_atoms, smiles)
        
        energy = self._compute_fake_energy(smiles, method)
        
        # Generate fake dipole moment (Debye)
        dipole = [
            random.gauss(0, 1.0),
            random.gauss(0, 1.0),
            random.gauss(0, 1.0),
        ]
        
        # Generate fake charges (Mulliken-like)
        charges = [random.gauss(0, 0.1) for _ in range(n_atoms)]
        # Normalize to zero total charge
        total_charge = sum(charges)
        charges = [c - total_charge / n_atoms for c in charges]
        
        return QMResult(
            energy=round(energy, 6),
            coordinates=coordinates,
            dipole=[round(d, 4) for d in dipole],
            charges=[round(c, 4) for c in charges],
            metadata={
                'method': method,
                'basis': basis,
                'smiles': smiles,
                'mock': True,
            }
        )

