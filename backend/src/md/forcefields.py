"""
Force Fields - Mock force field implementations

Provides interfaces and simple mock force fields for MD simulations.
"""

from typing import List, Protocol
import numpy as np
import logging

logger = logging.getLogger(__name__)


class ForceField(Protocol):
    """Protocol for force field implementations."""
    
    def compute_forces(
        self,
        coordinates: List[List[float]],
        smiles: str
    ) -> List[List[float]]:
        """
        Compute forces on atoms.
        
        Args:
            coordinates: [[x, y, z], ...] per atom
            smiles: SMILES string
        
        Returns:
            [[fx, fy, fz], ...] forces per atom
        """
        ...
    
    def compute_potential_energy(
        self,
        coordinates: List[List[float]],
        smiles: str
    ) -> float:
        """
        Compute potential energy.
        
        Args:
            coordinates: [[x, y, z], ...] per atom
            smiles: SMILES string
        
        Returns:
            Potential energy
        """
        ...


class SimpleForceField:
    """
    Simple mock force field.
    
    Uses harmonic potentials and Lennard-Jones-like interactions.
    """
    
    def __init__(self):
        self.k_bond = 100.0  # Bond spring constant
        self.k_angle = 10.0  # Angle spring constant
        self.epsilon = 0.1  # LJ epsilon
        self.sigma = 3.0  # LJ sigma
    
    def compute_forces(
        self,
        coordinates: List[List[float]],
        smiles: str
    ) -> List[List[float]]:
        """
        Compute forces using simple harmonic and LJ potentials.
        
        This is a mock implementation that generates realistic-looking forces.
        """
        n_atoms = len(coordinates)
        forces = [[0.0, 0.0, 0.0] for _ in range(n_atoms)]
        
        if n_atoms < 2:
            return forces
        
        # Compute pairwise interactions
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                # Distance vector
                r = np.array(coordinates[j]) - np.array(coordinates[i])
                r_mag = np.linalg.norm(r)
                
                if r_mag < 1e-10:
                    continue
                
                # Unit vector
                r_unit = r / r_mag
                
                # Simple force: harmonic + repulsion
                # F = -k * (r - r0) - repulsion
                r0 = 1.5  # Equilibrium distance
                force_mag = -self.k_bond * (r_mag - r0)
                
                # Add repulsion at short distances
                if r_mag < 1.0:
                    force_mag += 10.0 / (r_mag**2)
                
                # Force vector
                force = force_mag * r_unit
                
                # Newton's third law
                forces[i] = [f - fi for f, fi in zip(forces[i], force)]
                forces[j] = [f + fi for f, fi in zip(forces[j], force)]
        
        return forces
    
    def compute_potential_energy(
        self,
        coordinates: List[List[float]],
        smiles: str
    ) -> float:
        """
        Compute potential energy.
        
        Uses harmonic and LJ-like potentials.
        """
        n_atoms = len(coordinates)
        if n_atoms < 2:
            return 0.0
        
        energy = 0.0
        
        # Pairwise interactions
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                r = np.array(coordinates[j]) - np.array(coordinates[i])
                r_mag = np.linalg.norm(r)
                
                if r_mag < 1e-10:
                    continue
                
                # Harmonic potential
                r0 = 1.5
                energy += 0.5 * self.k_bond * (r_mag - r0)**2
                
                # LJ-like repulsion
                if r_mag < 1.0:
                    energy += self.epsilon * (self.sigma / r_mag)**12
        
        return float(energy)

