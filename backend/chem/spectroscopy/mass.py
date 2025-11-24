"""
Mass Spectrometry Fragmentation Engine
"""
from typing import Dict, List, Any

def predict_mass_spectrum(molecule: Dict[str, Any]) -> Dict[str, Any]:
    """
    Predict mass spectrum with fragmentation
    
    Returns:
        {
            "molecular_ion": M,
            "peaks": [
                {
                    "m/z": 43,
                    "intensity": 0.8,
                    "fragment_atoms": [atom_ids],
                    "fragment_type": "alpha_cleavage"
                }
            ]
        }
    """
    atoms = molecule.get("atoms", [])
    bonds = molecule.get("bonds", [])
    
    # Calculate molecular weight
    atomic_weights = {
        "H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999,
        "F": 18.998, "S": 32.065, "P": 30.974, "Cl": 35.453,
        "Br": 79.904, "I": 126.904
    }
    
    molecular_ion = sum(atomic_weights.get(a["element"], 0) for a in atoms)
    
    # Simple fragmentation rules
    peaks = []
    
    # Molecular ion
    peaks.append({
        "m/z": round(molecular_ion, 1),
        "intensity": 1.0,
        "fragment_atoms": [a["id"] for a in atoms],
        "fragment_type": "molecular_ion"
    })
    
    # Common fragments
    # M-15 (CH3 loss)
    peaks.append({
        "m/z": round(molecular_ion - 15.0, 1),
        "intensity": 0.3,
        "fragment_atoms": [],
        "fragment_type": "methyl_loss"
    })
    
    # M-17 (OH loss)
    peaks.append({
        "m/z": round(molecular_ion - 17.0, 1),
        "intensity": 0.2,
        "fragment_atoms": [],
        "fragment_type": "hydroxyl_loss"
    })
    
    # M-18 (H2O loss)
    peaks.append({
        "m/z": round(molecular_ion - 18.0, 1),
        "intensity": 0.4,
        "fragment_atoms": [],
        "fragment_type": "water_loss"
    })
    
    # Common small fragments
    common_fragments = [
        (15, "CH3"),
        (29, "C2H5"),
        (43, "C3H7"),
        (57, "C4H9"),
    ]
    
    for mz, name in common_fragments:
        peaks.append({
            "m/z": mz,
            "intensity": 0.5,
            "fragment_atoms": [],
            "fragment_type": name
        })
    
    # Sort by m/z
    peaks.sort(key=lambda x: x["m/z"], reverse=True)
    
    return {
        "molecular_ion": round(molecular_ion, 1),
        "peaks": peaks
    }

