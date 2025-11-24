"""
Retrosynthesis planning algorithm (graph search)
"""
from typing import Dict, List, Any, Optional, Tuple
import copy
from .utils import find_applicable_reactions, get_reaction_precursors
from chem.search.smarts import match_smarts

class PathwayNode:
    """Node in retrosynthesis pathway"""
    def __init__(self, molecule: Dict[str, Any], step: int, parent: Optional['PathwayNode'] = None):
        self.molecule = molecule
        self.step = step
        self.parent = parent
        self.reaction = None
        self.precursors = []
        self.cost = 0.0

def is_simple_molecule(molecule: Dict[str, Any]) -> bool:
    """
    Check if molecule is simple enough to be a starting material
    (e.g., common reagents like alcohols, acids)
    """
    atoms = molecule.get("atoms", [])
    bonds = molecule.get("bonds", [])
    
    # Simple heuristics: small molecules (< 5 atoms) or common patterns
    if len(atoms) <= 4:
        return True
    
    # Check for common starting materials (simplified)
    simple_patterns = ["[C][OH]", "[C]=O", "[C][NH2]"]
    for pattern in simple_patterns:
        is_match, _ = match_smarts(pattern, molecule)
        if is_match and len(atoms) <= 6:
            return True
    
    return False

def plan_retrosynthesis(
    target_molecule: Dict[str, Any],
    max_steps: int = 5,
    max_pathways: int = 10
) -> List[Dict[str, Any]]:
    """
    Plan retrosynthesis pathways using graph search (BFS)
    
    Returns:
        List of pathways, each with steps, precursors, and scores
    """
    pathways = []
    
    # BFS queue: (molecule, step, pathway_history)
    queue = [(target_molecule, 0, [])]
    visited = set()
    
    while queue and len(pathways) < max_pathways:
        current_mol, step, history = queue.pop(0)
        
        # Create molecule signature for visited check
        mol_sig = _molecule_signature(current_mol)
        if mol_sig in visited:
            continue
        visited.add(mol_sig)
        
        # Check if we've reached a simple starting material
        if is_simple_molecule(current_mol):
            pathway = {
                "steps": history + [{"molecule": current_mol, "step": step, "is_starting": True}],
                "total_steps": step,
                "target": target_molecule
            }
            pathways.append(pathway)
            continue
        
        # Stop if max steps reached
        if step >= max_steps:
            continue
        
        # Find applicable reactions
        applicable_reactions = find_applicable_reactions(current_mol)
        
        if not applicable_reactions:
            # No reactions found - treat as dead end or starting material
            if step > 0:  # Only if we've made progress
                pathway = {
                    "steps": history + [{"molecule": current_mol, "step": step, "is_starting": True}],
                    "total_steps": step,
                    "target": target_molecule
                }
                pathways.append(pathway)
            continue
        
        # For each applicable reaction, generate precursors
        for reaction in applicable_reactions[:3]:  # Limit to top 3 reactions per step
            precursors = get_reaction_precursors(reaction, current_mol)
            
            if not precursors:
                continue
            
            # Create new pathway step
            new_history = history + [{
                "molecule": current_mol,
                "step": step,
                "reaction": reaction,
                "precursors": precursors
            }]
            
            # Add each precursor to queue
            for precursor in precursors:
                if precursor.get("atoms"):
                    queue.append((precursor, step + 1, new_history))
    
    return pathways

def _molecule_signature(molecule: Dict[str, Any]) -> str:
    """Create a simple signature for molecule (for visited tracking)"""
    atoms = molecule.get("atoms", [])
    bonds = molecule.get("bonds", [])
    
    # Simple signature: element counts
    element_counts = {}
    for atom in atoms:
        elem = atom.get("element", "?")
        element_counts[elem] = element_counts.get(elem, 0) + 1
    
    sig = "_".join(f"{k}{v}" for k, v in sorted(element_counts.items()))
    return sig

def expand_pathway(pathway: Dict[str, Any], max_steps: int = 5) -> Dict[str, Any]:
    """
    Expand a pathway by continuing retrosynthesis from precursors
    """
    expanded_steps = []
    
    for step_data in pathway.get("steps", []):
        expanded_steps.append(step_data)
        
        precursors = step_data.get("precursors", [])
        if not precursors:
            continue
        
        # For each precursor, try to find further reactions
        for precursor in precursors:
            if is_simple_molecule(precursor):
                continue
            
            # Recursively plan from precursor (limited depth)
            sub_pathways = plan_retrosynthesis(precursor, max_steps=2, max_pathways=1)
            if sub_pathways:
                expanded_steps.extend(sub_pathways[0].get("steps", []))
    
    return {
        **pathway,
        "steps": expanded_steps,
        "total_steps": len(expanded_steps)
    }

