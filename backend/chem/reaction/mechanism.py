"""
Mechanism prediction - stepwise reaction intermediates
"""
from typing import Dict, List, Any
import copy
from .engine import predict_reaction_products, apply_reaction
from chem.energy.forcefield import calculate_total_energy

def predict_mechanism(
    reactants: List[Dict[str, Any]],
    reaction_type: str,
    max_steps: int = 5
) -> Dict[str, Any]:
    """
    Predict stepwise reaction mechanism
    
    Returns:
        {
            "steps": [
                {
                    "step": int,
                    "intermediate": molecule_json,
                    "energy": float,
                    "description": str
                }
            ],
            "activation_energy": float,
            "reaction_energy": float
        }
    """
    steps = []
    
    # Initial state
    if len(reactants) > 0:
        initial = reactants[0]
        initial_energy = calculate_total_energy(initial)["total_energy"]
        steps.append({
            "step": 0,
            "intermediate": copy.deepcopy(initial),
            "energy": initial_energy,
            "description": "Reactant"
        })
    
    # Simulate reaction steps
    current_molecule = copy.deepcopy(reactants[0]) if reactants else None
    if not current_molecule:
        return {"steps": [], "activation_energy": 0.0, "reaction_energy": 0.0}
    
    for step in range(1, max_steps + 1):
        # Try to apply reaction
        products = apply_reaction(current_molecule, reaction_type)
        
        if not products:
            break
        
        # Take first product
        product = products[0]
        product_energy = calculate_total_energy(product)["total_energy"]
        
        steps.append({
            "step": step,
            "intermediate": copy.deepcopy(product),
            "energy": product_energy,
            "description": f"Intermediate {step}"
        })
        
        current_molecule = product
    
    # Calculate activation and reaction energies
    if len(steps) >= 2:
        energies = [s["energy"] for s in steps]
        max_energy = max(energies)
        initial_energy = steps[0]["energy"]
        final_energy = steps[-1]["energy"]
        
        activation_energy = max_energy - initial_energy
        reaction_energy = final_energy - initial_energy
    else:
        activation_energy = 0.0
        reaction_energy = 0.0
    
    return {
        "steps": steps,
        "activation_energy": activation_energy,
        "reaction_energy": reaction_energy
    }

