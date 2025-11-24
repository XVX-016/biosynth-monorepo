"""
Geometry optimization using gradient descent
"""
from typing import Dict, List, Any
import copy
import math
from .forcefield import calculate_total_energy, distance

def calculate_gradient(molecule: Dict[str, Any], atom_id: str) -> List[float]:
    """
    Calculate energy gradient for a specific atom
    Returns [gx, gy, gz]
    """
    delta = 0.001  # Small displacement for numerical gradient
    
    atoms = {a["id"]: a for a in molecule["atoms"]}
    atom = atoms[atom_id]
    original_pos = atom["position"][:]
    
    gradient = [0.0, 0.0, 0.0]
    
    for axis in range(3):
        # Forward difference
        atom["position"][axis] = original_pos[axis] + delta
        energy_plus = calculate_total_energy(molecule)["total_energy"]
        
        # Backward difference
        atom["position"][axis] = original_pos[axis] - delta
        energy_minus = calculate_total_energy(molecule)["total_energy"]
        
        # Central difference
        gradient[axis] = (energy_plus - energy_minus) / (2 * delta)
        
        # Restore position
        atom["position"][axis] = original_pos[axis]
    
    return gradient

def optimize_geometry(
    molecule: Dict[str, Any],
    max_iterations: int = 100,
    convergence_threshold: float = 0.001,
    step_size: float = 0.01
) -> Dict[str, Any]:
    """
    Optimize molecular geometry using steepest descent
    
    Returns:
        {
            "optimized_molecule": molecule_json,
            "final_energy": float,
            "iterations": int,
            "energy_history": [float]
        }
    """
    mol = copy.deepcopy(molecule)
    atoms = {a["id"]: a for a in mol["atoms"]}
    energy_history = []
    
    initial_energy = calculate_total_energy(mol)["total_energy"]
    energy_history.append(initial_energy)
    
    for iteration in range(max_iterations):
        # Calculate gradients for all atoms
        gradients = {}
        for atom_id in atoms.keys():
            gradients[atom_id] = calculate_gradient(mol, atom_id)
        
        # Update positions
        max_gradient = 0.0
        for atom_id, grad in gradients.items():
            atom = atoms[atom_id]
            # Steepest descent: x_new = x_old - step_size * gradient
            for axis in range(3):
                atom["position"][axis] -= step_size * grad[axis]
            
            grad_mag = math.sqrt(sum(g*g for g in grad))
            max_gradient = max(max_gradient, grad_mag)
        
        # Check convergence
        current_energy = calculate_total_energy(mol)["total_energy"]
        energy_history.append(current_energy)
        
        if max_gradient < convergence_threshold:
            break
        
        # Adaptive step size
        if len(energy_history) > 1:
            if current_energy > energy_history[-2]:
                step_size *= 0.5  # Reduce step if energy increased
            else:
                step_size *= 1.05  # Increase step if energy decreased
        
        step_size = max(0.001, min(0.1, step_size))  # Clamp step size
    
    final_energy = calculate_total_energy(mol)["total_energy"]
    
    return {
        "optimized_molecule": mol,
        "final_energy": final_energy,
        "initial_energy": initial_energy,
        "iterations": iteration + 1,
        "energy_history": energy_history
    }

