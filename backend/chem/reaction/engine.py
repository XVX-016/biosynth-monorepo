"""
Reaction simulation engine
"""
from typing import Dict, List, Any, Tuple, Optional
import copy
from chem.search.smarts import match_smarts

# Common reaction patterns (SMARTS)
REACTION_PATTERNS = {
    "substitution": {
        "pattern": "[C:1][Cl:2]",
        "product_pattern": "[C:1][OH:2]",
        "type": "nucleophilic_substitution"
    },
    "addition": {
        "pattern": "[C:1]=[C:2]",
        "product_pattern": "[C:1][C:2]",
        "type": "addition"
    },
    "elimination": {
        "pattern": "[C:1][C:2][OH:3]",
        "product_pattern": "[C:1]=[C:2]",
        "type": "elimination"
    }
}

def find_reaction_sites(molecule: Dict[str, Any], pattern: str) -> List[Dict[str, Any]]:
    """
    Find atoms that match a reaction pattern
    """
    from chem.search.smarts import match_smarts
    
    is_match, matched_atoms = match_smarts(pattern, molecule)
    if is_match and matched_atoms:
        return matched_atoms
    return []

def apply_reaction(
    molecule: Dict[str, Any],
    reaction_type: str,
    conditions: Optional[Dict[str, Any]] = None
) -> List[Dict[str, Any]]:
    """
    Apply a reaction to a molecule
    
    Returns list of possible products
    """
    if reaction_type not in REACTION_PATTERNS:
        return []
    
    pattern_info = REACTION_PATTERNS[reaction_type]
    sites = find_reaction_sites(molecule, pattern_info["pattern"])
    
    if not sites:
        return []
    
    products = []
    for site in sites:
        product = simulate_reaction_step(molecule, site, pattern_info)
        if product:
            products.append(product)
    
    return products

def simulate_reaction_step(
    molecule: Dict[str, Any],
    reaction_site: Dict[str, Any],
    pattern_info: Dict[str, Any]
) -> Optional[Dict[str, Any]]:
    """
    Simulate a single reaction step
    """
    product = copy.deepcopy(molecule)
    
    # Simple bond modification based on pattern
    # In a real implementation, this would use SMARTS transformation rules
    atoms = {a["id"]: a for a in product["atoms"]}
    bonds = product["bonds"]
    
    # Example: substitution reaction - break one bond, form another
    if pattern_info["type"] == "nucleophilic_substitution":
        # Find the bond to break (C-Cl)
        for bond in bonds:
            a1 = atoms.get(bond["atom1"])
            a2 = atoms.get(bond["atom2"])
            if a1 and a2:
                if (a1["element"] == "C" and a2["element"] == "Cl") or \
                   (a1["element"] == "Cl" and a2["element"] == "C"):
                    # Remove Cl atom and bond
                    product["atoms"] = [a for a in product["atoms"] if a["id"] != a2["id"]]
                    product["bonds"] = [b for b in product["bonds"] if b["id"] != bond["id"]]
                    
                    # Add OH group
                    import uuid
                    o_id = str(uuid.uuid4())[:8]
                    h_id = str(uuid.uuid4())[:8]
                    
                    c_pos = a1["position"]
                    o_pos = [c_pos[0] + 1.4, c_pos[1], c_pos[2]]
                    h_pos = [o_pos[0] + 0.96, o_pos[1], o_pos[2]]
                    
                    product["atoms"].extend([
                        {"id": o_id, "element": "O", "position": o_pos},
                        {"id": h_id, "element": "H", "position": h_pos}
                    ])
                    import uuid
                    product["bonds"].extend([
                        {"id": str(uuid.uuid4())[:8], "atom1": a1["id"], "atom2": o_id, "order": 1},
                        {"id": str(uuid.uuid4())[:8], "atom1": o_id, "atom2": h_id, "order": 1}
                    ])
                    
                    return product
    
    return None

def predict_reaction_products(
    reactants: List[Dict[str, Any]],
    reaction_type: str,
    conditions: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    Predict products from reactants
    
    Returns:
        {
            "products": [molecule_json],
            "reaction_type": str,
            "feasible": bool
        }
    """
    if len(reactants) == 0:
        return {"products": [], "reaction_type": reaction_type, "feasible": False}
    
    # For now, apply reaction to first reactant
    main_reactant = reactants[0]
    products = apply_reaction(main_reactant, reaction_type, conditions)
    
    return {
        "products": products,
        "reaction_type": reaction_type,
        "feasible": len(products) > 0
    }

