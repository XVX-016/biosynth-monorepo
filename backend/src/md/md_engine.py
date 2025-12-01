"""
MD Engine - Mock molecular dynamics engine

Provides mock implementations of MD simulations.
All results are fake but follow realistic patterns.
"""

from typing import List, Optional, Dict, Any
import numpy as np
import random
import logging

from .md_interfaces import MDEngineProtocol, MDResult, TrajectoryFrame
from .forcefields import ForceField, SimpleForceField
from .integrators import Integrator, VerletIntegrator

logger = logging.getLogger(__name__)


class MDEngine:
    """
    Mock molecular dynamics engine.
    
    Returns fake but realistic-looking trajectories.
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
        return sum(1 for c in smiles if c.isupper() and c.isalpha())
    
    def _generate_initial_coordinates(
        self,
        n_atoms: int,
        smiles: str
    ) -> List[List[float]]:
        """Generate initial 3D coordinates."""
        coords = []
        radius = 2.0
        
        for i in range(n_atoms):
            theta = 2 * np.pi * i / n_atoms if n_atoms > 0 else 0
            phi = np.pi * random.random()
            
            x = radius * np.sin(phi) * np.cos(theta)
            y = radius * np.sin(phi) * np.sin(theta)
            z = radius * np.cos(phi)
            
            coords.append([float(x), float(y), float(z)])
        
        return coords
    
    def _compute_kinetic_energy(
        self,
        velocities: List[List[float]],
        masses: List[float]
    ) -> float:
        """Compute kinetic energy from velocities."""
        ke = 0.0
        for v, m in zip(velocities, masses):
            v_mag = np.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
            ke += 0.5 * m * v_mag**2
        return ke
    
    def _compute_temperature(
        self,
        kinetic_energy: float,
        n_atoms: int
    ) -> float:
        """Compute temperature from kinetic energy."""
        # T = 2 * KE / (3 * N * kB)
        # kB = 0.001987 kcal/(mol*K) ≈ 0.002
        kB = 0.002
        if n_atoms == 0:
            return 0.0
        return 2 * kinetic_energy / (3 * n_atoms * kB)
    
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
        Run MD simulation (mock).
        
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
        self._set_seed(smiles)
        
        n_atoms = self._estimate_atom_count(smiles)
        
        # Generate initial coordinates
        if coordinates is None:
            coords = self._generate_initial_coordinates(n_atoms, smiles)
        else:
            coords = [list(c) for c in coordinates]
        
        # Initialize force field and integrator
        ff = SimpleForceField()
        integrator = VerletIntegrator(timestep)
        
        # Initialize velocities (random, scaled to target temperature)
        masses = [12.0] * n_atoms  # Assume carbon-like atoms
        velocities = []
        for _ in range(n_atoms):
            # Maxwell-Boltzmann distribution
            v = [
                random.gauss(0, np.sqrt(temperature / 100.0)),
                random.gauss(0, np.sqrt(temperature / 100.0)),
                random.gauss(0, np.sqrt(temperature / 100.0)),
            ]
            velocities.append(v)
        
        # Generate trajectory
        trajectory = []
        current_coords = coords
        current_velocities = velocities
        
        for step in range(steps):
            # Compute forces (mock)
            forces = ff.compute_forces(current_coords, smiles)
            
            # Integrate (mock Verlet)
            current_coords, current_velocities = integrator.step(
                current_coords,
                current_velocities,
                forces,
                masses
            )
            
            # Compute energy (mock)
            ke = self._compute_kinetic_energy(current_velocities, masses)
            pe = ff.compute_potential_energy(current_coords, smiles)
            total_energy = ke + pe
            
            # Compute temperature
            temp = self._compute_temperature(ke, n_atoms)
            
            # Store frame (every 10 steps to reduce size)
            if step % 10 == 0 or step == steps - 1:
                trajectory.append(TrajectoryFrame(
                    step=step,
                    time=step * timestep,
                    coordinates=[list(c) for c in current_coords],
                    velocities=[list(v) for v in current_velocities],
                    energy=round(total_energy, 6),
                    temperature=round(temp, 2),
                ))
        
        return MDResult(
            trajectory=trajectory,
            final_coordinates=[list(c) for c in current_coords],
            total_energy=round(total_energy, 6),
            metadata={
                'smiles': smiles,
                'steps': steps,
                'timestep': timestep,
                'temperature': temperature,
                'forcefield': forcefield,
                'mock': True,
            }
        )

