"""
Reaction database access and utilities for retrosynthesis
"""
from typing import Dict, List, Any, Optional
from chem.search.smarts import match_smarts

# Simplified reaction database (in production, this would be a real DB)
REACTION_DATABASE = [
    {
        "id": "rxn_001",
        "name": "Nucleophilic Substitution",
        "product_pattern": "[C:1][Cl:2]",
        "reactant_patterns": [
            {"pattern": "[C:1][OH:2]", "reagent": "NaOH"},
            {"pattern": "[C:1][NH2:2]", "reagent": "NH3"}
        ],
        "conditions": {"temperature": 25, "solvent": "water"},
        "yield": 0.85,
        "cost": 10.0
    },
    {
        "id": "rxn_002",
        "name": "Esterification",
        "product_pattern": "[C:1][O:2][C:3]=O",
        "reactant_patterns": [
            {"pattern": "[C:1][OH:2]", "reagent": "Acid + Alcohol"},
            {"pattern": "[C:3][OH:4]", "reagent": "Acid + Alcohol"}
        ],
        "conditions": {"temperature": 80, "solvent": "organic"},
        "yield": 0.75,
        "cost": 15.0
    },
    {
        "id": "rxn_003",
        "name": "Reduction",
        "product_pattern": "[C:1]=[O:2]",
        "reactant_patterns": [
            {"pattern": "[C:1][OH:2]", "reagent": "LiAlH4"},
            {"pattern": "[C:1][H:2]", "reagent": "NaBH4"}
        ],
        "conditions": {"temperature": 0, "solvent": "THF"},
        "yield": 0.90,
        "cost": 25.0
    },
    {
        "id": "rxn_004",
        "name": "Oxidation",
        "product_pattern": "[C:1][OH:2]",
        "reactant_patterns": [
            {"pattern": "[C:1]=[O:2]", "reagent": "KMnO4"},
            {"pattern": "[C:1][H:2]", "reagent": "CrO3"}
        ],
        "conditions": {"temperature": 50, "solvent": "water"},
        "yield": 0.80,
        "cost": 20.0
    }
]

def find_applicable_reactions(target_molecule: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Find reactions that can produce the target molecule
    
    Returns list of applicable reactions with matched patterns
    """
    applicable = []
    
    for reaction in REACTION_DATABASE:
        product_pattern = reaction.get("product_pattern", "")
        if not product_pattern:
            continue
        
        # Check if target matches product pattern
        is_match, matched_atoms = match_smarts(product_pattern, target_molecule)
        if is_match:
            applicable.append({
                **reaction,
                "matched_atoms": matched_atoms
            })
    
    return applicable

def get_reaction_precursors(reaction: Dict[str, Any], target_molecule: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Generate precursor molecules for a given reaction
    
    Returns list of precursor molecule structures
    """
    precursors = []
    reactant_patterns = reaction.get("reactant_patterns", [])
    
    for reactant_info in reactant_patterns:
        pattern = reactant_info.get("pattern", "")
        reagent = reactant_info.get("reagent", "")
        
        # In a real implementation, this would use SMARTS transformation
        # For now, return simplified precursor structures
        # This is a placeholder - real retrosynthesis would apply inverse SMARTS
        
        # Simplified: create a basic precursor structure
        precursor = {
            "atoms": [],
            "bonds": [],
            "reagent": reagent,
            "reaction_id": reaction["id"]
        }
        
        # Try to extract atoms from pattern (simplified)
        # In production, use proper SMARTS inverse transformation
        if "[C:1]" in pattern and "[OH:2]" in pattern:
            # Example: precursor for substitution might be alcohol
            precursor["atoms"] = [
                {"id": "c1", "element": "C", "position": [0, 0, 0]},
                {"id": "o1", "element": "O", "position": [1.4, 0, 0]},
                {"id": "h1", "element": "H", "position": [2.0, 0, 0]}
            ]
            precursor["bonds"] = [
                {"id": "b1", "atom1": "c1", "atom2": "o1", "order": 1},
                {"id": "b2", "atom1": "o1", "atom2": "h1", "order": 1}
            ]
        
        if precursor["atoms"]:
            precursors.append(precursor)
    
    return precursors

def get_reagent_availability(reagent: str) -> Dict[str, Any]:
    """
    Check reagent availability and cost
    
    Returns availability info
    """
    # Simplified availability database
    availability_db = {
        "NaOH": {"available": True, "cost": 5.0, "supplier": "Sigma"},
        "NH3": {"available": True, "cost": 8.0, "supplier": "Sigma"},
        "LiAlH4": {"available": True, "cost": 30.0, "supplier": "TCI"},
        "NaBH4": {"available": True, "cost": 15.0, "supplier": "Sigma"},
        "KMnO4": {"available": True, "cost": 12.0, "supplier": "Sigma"},
        "CrO3": {"available": True, "cost": 18.0, "supplier": "Sigma"},
        "Acid + Alcohol": {"available": True, "cost": 10.0, "supplier": "Generic"}
    }
    
    return availability_db.get(reagent, {"available": False, "cost": 100.0, "supplier": "Unknown"})

