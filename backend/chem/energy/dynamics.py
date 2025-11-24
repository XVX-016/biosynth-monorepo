"""
Molecular dynamics simulation
"""
from typing import Dict, List, Any
import copy
import math
from .forcefield import calculate_total_energy, distance

def verlet_integration(
    molecule: Dict[str, Any],
    velocities: Dict[str, List[float]],
    forces: Dict[str, List[float]],
    dt: float
) -> Dict[str, Any]:
    """
    Verlet integration step
    Returns updated positions and velocities
    """
    atoms = {a["id"]: a for a in molecule["atoms"]}
    new_velocities = {}
    
    # Atomic masses (simplified)
    masses = {
        "H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999,
        "F": 18.998, "S": 32.065, "P": 30.974
    }
    
    for atom_id, atom in atoms.items():
        mass = masses.get(atom["element"], 12.0)
        force = forces.get(atom_id, [0.0, 0.0, 0.0])
        vel = velocities.get(atom_id, [0.0, 0.0, 0.0])
        
        # Update velocity: v = v + (F/m) * dt
        new_vel = [
            vel[0] + (force[0] / mass) * dt,
            vel[1] + (force[1] / mass) * dt,
            vel[2] + (force[2] / mass) * dt
        ]
        new_velocities[atom_id] = new_vel
        
        # Update position: x = x + v * dt
        atom["position"][0] += new_vel[0] * dt
        atom["position"][1] += new_vel[1] * dt
        atom["position"][2] += new_vel[2] * dt
    
    return new_velocities

def calculate_forces(molecule: Dict[str, Any]) -> Dict[str, List[float]]:
    """
    Calculate forces on all atoms (negative gradient of energy)
    """
    delta = 0.001
    atoms = {a["id"]: a for a in molecule["atoms"]}
    forces = {}
    
    for atom_id in atoms.keys():
        atom = atoms[atom_id]
        original_pos = atom["position"][:]
        force = [0.0, 0.0, 0.0]
        
        for axis in range(3):
            atom["position"][axis] = original_pos[axis] + delta
            energy_plus = calculate_total_energy(molecule)["total_energy"]
            
            atom["position"][axis] = original_pos[axis] - delta
            energy_minus = calculate_total_energy(molecule)["total_energy"]
            
            # Force = -dE/dx
            force[axis] = -(energy_plus - energy_minus) / (2 * delta)
            atom["position"][axis] = original_pos[axis]
        
        forces[atom_id] = force
    
    return forces

def scale_velocities(velocities: Dict[str, List[float]], target_temp: float) -> Dict[str, List[float]]:
    """
    Scale velocities to match target temperature (NVT ensemble)
    """
    masses = {
        "H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999,
        "F": 18.998, "S": 32.065, "P": 30.974
    }
    
    # Calculate current kinetic energy
    ke = 0.0
    for atom_id, vel in velocities.items():
        # Find atom to get mass
        # Simplified: assume average mass
        mass = 12.0  # Default
        v_mag_sq = vel[0]**2 + vel[1]**2 + vel[2]**2
        ke += 0.5 * mass * v_mag_sq
    
    # Target kinetic energy: KE = (3/2) * N * k_B * T
    # k_B = 0.001987 kcal/(mol·K), T in Kelvin
    kb = 0.001987
    n_atoms = len(velocities)
    target_ke = 1.5 * n_atoms * kb * target_temp
    
    if ke > 0:
        scale = math.sqrt(target_ke / ke)
        scaled_velocities = {}
        for atom_id, vel in velocities.items():
            scaled_velocities[atom_id] = [v * scale for v in vel]
        return scaled_velocities
    
    return velocities

def simulate_dynamics(
    molecule: Dict[str, Any],
    temperature: float = 300.0,
    steps: int = 100,
    dt: float = 0.001  # Time step in ps
) -> Dict[str, Any]:
    """
    Run molecular dynamics simulation
    
    Returns:
        {
            "trajectory": [ { "atoms": [...], "time": float } ],
            "energies": [ { "time": float, "energy": float } ],
            "final_molecule": molecule_json
        }
    """
    mol = copy.deepcopy(molecule)
    atoms = {a["id"]: a for a in mol["atoms"]}
    
    # Initialize velocities (random, then scale to temperature)
    import random
    velocities = {}
    for atom_id in atoms.keys():
        velocities[atom_id] = [
            random.gauss(0, 1),
            random.gauss(0, 1),
            random.gauss(0, 1)
        ]
    
    velocities = scale_velocities(velocities, temperature)
    
    trajectory = []
    energies = []
    
    # Initial snapshot
    trajectory.append({
        "atoms": copy.deepcopy(mol["atoms"]),
        "time": 0.0
    })
    initial_energy = calculate_total_energy(mol)["total_energy"]
    energies.append({
        "time": 0.0,
        "energy": initial_energy
    })
    
    # MD loop
    for step in range(steps):
        # Calculate forces
        forces = calculate_forces(mol)
        
        # Verlet integration
        velocities = verlet_integration(mol, velocities, forces, dt)
        
        # Periodic velocity scaling (NVT)
        if step % 10 == 0:
            velocities = scale_velocities(velocities, temperature)
        
        # Save snapshot
        if step % 10 == 0:  # Save every 10 steps
            time = (step + 1) * dt
            trajectory.append({
                "atoms": copy.deepcopy(mol["atoms"]),
                "time": time
            })
            energy = calculate_total_energy(mol)["total_energy"]
            energies.append({
                "time": time,
                "energy": energy
            })
    
    return {
        "trajectory": trajectory,
        "energies": energies,
        "final_molecule": mol
    }

