"""
IR Spectrum Prediction Engine
"""
from typing import Dict, List, Any
from .utils import detect_functional_groups

# IR frequency ranges (cm^-1)
IR_RANGES = {
    "O-H": (3200, 3600, "broad"),
    "N-H": (3300, 3500, "medium"),
    "C-H_aromatic": (3000, 3100, "weak"),
    "C-H_aliphatic": (2850, 3000, "medium"),
    "C=O_ketone": (1700, 1725, "strong"),
    "C=O_aldehyde": (1720, 1740, "strong"),
    "C=O_carboxylic_acid": (1710, 1720, "strong"),
    "C=O_ester": (1735, 1750, "strong"),
    "C=C": (1600, 1680, "medium"),
    "C=C_aromatic": (1450, 1600, "medium"),
    "C-O": (1000, 1300, "strong"),
    "C-N": (1000, 1350, "medium"),
}

def predict_ir_spectrum(molecule: Dict[str, Any]) -> Dict[str, Any]:
    """
    Predict IR spectrum for molecule
    
    Returns:
        {
            "peaks": [
                {
                    "wavenumber": 1715,
                    "intensity": "strong",
                    "group": "C=O",
                    "atoms": [atom_ids]
                }
            ]
        }
    """
    groups = detect_functional_groups(molecule)
    peaks = []
    
    for group in groups:
        group_type = group["type"]
        subtype = group.get("subtype", "")
        
        if group_type == "carbonyl":
            key = f"C=O_{subtype}" if subtype in ["ketone", "aldehyde", "carboxylic_acid", "ester"] else "C=O_ketone"
            if key in IR_RANGES:
                low, high, intensity = IR_RANGES[key]
                wavenumber = (low + high) / 2
                peaks.append({
                    "wavenumber": round(wavenumber, 1),
                    "intensity": intensity,
                    "group": f"C=O ({subtype})",
                    "atoms": group["atoms"]
                })
        
        elif group_type == "hydroxyl":
            low, high, intensity = IR_RANGES["O-H"]
            wavenumber = (low + high) / 2
            peaks.append({
                "wavenumber": round(wavenumber, 1),
                "intensity": intensity,
                "group": "O-H",
                "atoms": group["atoms"]
            })
        
        elif group_type == "alkene":
            low, high, intensity = IR_RANGES["C=C"]
            wavenumber = (low + high) / 2
            peaks.append({
                "wavenumber": round(wavenumber, 1),
                "intensity": intensity,
                "group": "C=C",
                "atoms": group["atoms"]
            })
        
        elif group_type == "aromatic":
            low, high, intensity = IR_RANGES["C=C_aromatic"]
            wavenumber = (low + high) / 2
            peaks.append({
                "wavenumber": round(wavenumber, 1),
                "intensity": intensity,
                "group": "Aromatic C=C",
                "atoms": group["atoms"]
            })
    
    # Add C-H stretches
    atoms = molecule.get("atoms", [])
    has_aromatic = any(g["type"] == "aromatic" for g in groups)
    for atom in atoms:
        if atom["element"] == "H":
            if has_aromatic:
                low, high, intensity = IR_RANGES["C-H_aromatic"]
            else:
                low, high, intensity = IR_RANGES["C-H_aliphatic"]
            wavenumber = (low + high) / 2
            peaks.append({
                "wavenumber": round(wavenumber, 1),
                "intensity": intensity,
                "group": "C-H",
                "atoms": [atom["id"]]
            })
    
    # Sort by wavenumber
    peaks.sort(key=lambda x: x["wavenumber"], reverse=True)
    
    return {"peaks": peaks}

