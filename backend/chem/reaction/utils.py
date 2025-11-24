"""
Reaction thermodynamics utilities
"""
from typing import Dict, List, Any
from chem.energy.forcefield import calculate_total_energy

def calculate_reaction_energy(
    reactants: List[Dict[str, Any]],
    products: List[Dict[str, Any]]
) -> Dict[str, float]:
    """
    Calculate reaction thermodynamics
    
    Returns:
        {
            "delta_E": float,
            "delta_H": float,  # Approximate
            "delta_G": float  # Approximate
        }
    """
    reactant_energy = sum(calculate_total_energy(r)["total_energy"] for r in reactants)
    product_energy = sum(calculate_total_energy(p)["total_energy"] for p in products)
    
    delta_E = product_energy - reactant_energy
    
    # Approximate ΔH ≈ ΔE (for gas phase, no significant volume change)
    delta_H = delta_E
    
    # Approximate ΔG ≈ ΔH (neglecting entropy for now)
    delta_G = delta_H
    
    return {
        "delta_E": delta_E,
        "delta_H": delta_H,
        "delta_G": delta_G
    }

def find_activation_barrier(mechanism: Dict[str, Any]) -> float:
    """
    Find activation energy barrier from mechanism
    """
    steps = mechanism.get("steps", [])
    if len(steps) < 2:
        return 0.0
    
    energies = [s["energy"] for s in steps]
    initial_energy = energies[0]
    max_energy = max(energies)
    
    return max_energy - initial_energy

