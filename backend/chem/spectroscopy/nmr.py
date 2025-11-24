"""
NMR Spectrum Prediction Engine (1H and 13C)
"""
from typing import Dict, List, Any
from .utils import detect_functional_groups

# Chemical shift ranges (ppm)
H_SHIFTS = {
    "alkane": (0.8, 1.5),
    "alkene": (4.5, 6.5),
    "aromatic": (6.5, 8.5),
    "aldehyde": (9.5, 10.5),
    "carboxylic_acid": (11.0, 12.0),
    "alcohol": (1.0, 5.0),
    "amine": (1.0, 5.0),
}

C_SHIFTS = {
    "alkane": (10, 60),
    "alkene": (100, 150),
    "aromatic": (110, 150),
    "carbonyl_ketone": (200, 220),
    "carbonyl_aldehyde": (190, 205),
    "carbonyl_acid": (165, 185),
    "carbonyl_ester": (165, 175),
}

def predict_nmr_spectrum(molecule: Dict[str, Any]) -> Dict[str, Any]:
    """
    Predict 1H and 13C NMR spectra
    
    Returns:
        {
            "protons": [
                {
                    "shift": 7.2,
                    "multiplicity": "doublet",
                    "integration": 2,
                    "atom": atom_id
                }
            ],
            "carbons": [
                {
                    "shift": 120.5,
                    "atom": atom_id
                }
            ]
        }
    """
    atoms = molecule.get("atoms", [])
    bonds = molecule.get("bonds", [])
    groups = detect_functional_groups(molecule)
    
    # Build adjacency
    bonds_by_atom: Dict[str, List[Dict]] = {a["id"]: [] for a in atoms}
    for bond in bonds:
        bonds_by_atom[bond["atom1"]].append(bond)
        bonds_by_atom[bond["atom2"]].append(bond)
    
    protons = []
    carbons = []
    
    has_aromatic = any(g["type"] == "aromatic" for g in groups)
    
    for atom in atoms:
        element = atom["element"]
        atom_id = atom["id"]
        neighbors = bonds_by_atom.get(atom_id, [])
        
        if element == "H":
            # Determine chemical environment
            if has_aromatic:
                low, high = H_SHIFTS["aromatic"]
                multiplicity = "singlet"  # Simplified
            else:
                # Check if attached to carbonyl
                attached_to_carbonyl = False
                for bond in neighbors:
                    other_id = bond["atom1"] if bond["atom1"] != atom_id else bond["atom2"]
                    other_atom = next((a for a in atoms if a["id"] == other_id), None)
                    if other_atom:
                        for group in groups:
                            if group["type"] == "carbonyl" and other_id in group["atoms"]:
                                attached_to_carbonyl = True
                                break
                
                if attached_to_carbonyl:
                    low, high = H_SHIFTS["aldehyde"]
                    multiplicity = "singlet"
                else:
                    low, high = H_SHIFTS["alkane"]
                    multiplicity = "singlet"
            
            shift = (low + high) / 2
            protons.append({
                "shift": round(shift, 2),
                "multiplicity": multiplicity,
                "integration": 1,
                "atom": atom_id
            })
        
        elif element == "C":
            # Determine chemical environment
            if has_aromatic:
                low, high = C_SHIFTS["aromatic"]
            else:
                # Check if carbonyl
                is_carbonyl = False
                carbonyl_type = None
                for group in groups:
                    if group["type"] == "carbonyl" and atom_id in group["atoms"]:
                        is_carbonyl = True
                        carbonyl_type = group.get("subtype", "ketone")
                        break
                
                if is_carbonyl:
                    key = f"carbonyl_{carbonyl_type}" if carbonyl_type in ["ketone", "aldehyde", "carboxylic_acid", "ester"] else "carbonyl_ketone"
                    if key in C_SHIFTS:
                        low, high = C_SHIFTS[key]
                    else:
                        low, high = C_SHIFTS["carbonyl_ketone"]
                else:
                    # Check for double bond
                    has_double_bond = any(b.get("order", 1) == 2 for b in neighbors)
                    if has_double_bond:
                        low, high = C_SHIFTS["alkene"]
                    else:
                        low, high = C_SHIFTS["alkane"]
            
            shift = (low + high) / 2
            carbons.append({
                "shift": round(shift, 1),
                "atom": atom_id
            })
    
    # Sort by shift
    protons.sort(key=lambda x: x["shift"], reverse=True)
    carbons.sort(key=lambda x: x["shift"], reverse=True)
    
    return {
        "protons": protons,
        "carbons": carbons
    }

