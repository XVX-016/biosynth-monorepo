"""
Functional group detection utilities
"""
from typing import Dict, List, Any, Set, Tuple

def detect_functional_groups(molecule: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Detect functional groups in molecule
    
    Returns:
        [
            {
                "type": "carbonyl",
                "subtype": "ketone",
                "atoms": [atom_ids],
                "bonds": [bond_ids]
            },
            ...
        ]
    """
    atoms = molecule.get("atoms", [])
    bonds = molecule.get("bonds", [])
    
    # Build adjacency
    atom_map = {a["id"]: a for a in atoms}
    bonds_by_atom: Dict[str, List[Dict]] = {a["id"]: [] for a in atoms}
    for bond in bonds:
        bonds_by_atom[bond["atom1"]].append(bond)
        bonds_by_atom[bond["atom2"]].append(bond)
    
    groups = []
    
    # Detect C=O (carbonyl)
    for bond in bonds:
        if bond.get("order", 1) == 2:
            a1 = atom_map.get(bond["atom1"])
            a2 = atom_map.get(bond["atom2"])
            if a1 and a2:
                if (a1["element"] == "C" and a2["element"] == "O") or \
                   (a1["element"] == "O" and a2["element"] == "C"):
                    carbon = a1 if a1["element"] == "C" else a2
                    oxygen = a2 if a1["element"] == "C" else a1
                    
                    # Determine subtype
                    carbon_bonds = bonds_by_atom.get(carbon["id"], [])
                    has_h = any(atom_map.get(b["atom1"] if b["atom1"] != carbon["id"] else b["atom2"], {}).get("element") == "H" 
                               for b in carbon_bonds if b["atom1"] == carbon["id"] or b["atom2"] == carbon["id"])
                    has_oh = any(atom_map.get(b["atom1"] if b["atom1"] != carbon["id"] else b["atom2"], {}).get("element") == "O"
                               for b in carbon_bonds if b["atom1"] == carbon["id"] or b["atom2"] == carbon["id"])
                    
                    subtype = "aldehyde" if has_h else "carboxylic_acid" if has_oh else "ketone"
                    
                    groups.append({
                        "type": "carbonyl",
                        "subtype": subtype,
                        "atoms": [carbon["id"], oxygen["id"]],
                        "bonds": [bond["id"]]
                    })
    
    # Detect OH groups
    for atom in atoms:
        if atom["element"] == "O":
            o_bonds = bonds_by_atom.get(atom["id"], [])
            h_bonds = [b for b in o_bonds if atom_map.get(b["atom1"] if b["atom1"] != atom["id"] else b["atom2"], {}).get("element") == "H"]
            if h_bonds:
                groups.append({
                    "type": "hydroxyl",
                    "subtype": "alcohol",
                    "atoms": [atom["id"], h_bonds[0]["atom1"] if h_bonds[0]["atom1"] != atom["id"] else h_bonds[0]["atom2"]],
                    "bonds": [h_bonds[0]["id"]]
                })
    
    # Detect C=C (alkene)
    for bond in bonds:
        if bond.get("order", 1) == 2:
            a1 = atom_map.get(bond["atom1"])
            a2 = atom_map.get(bond["atom2"])
            if a1 and a2 and a1["element"] == "C" and a2["element"] == "C":
                groups.append({
                    "type": "alkene",
                    "subtype": "double_bond",
                    "atoms": [a1["id"], a2["id"]],
                    "bonds": [bond["id"]]
                })
    
    # Detect aromatic rings (simplified - 6-membered rings)
    aromatic_atoms = _detect_aromatic_rings(atoms, bonds, atom_map, bonds_by_atom)
    if aromatic_atoms:
        groups.append({
            "type": "aromatic",
            "subtype": "ring",
            "atoms": aromatic_atoms,
            "bonds": []
        })
    
    return groups

def _detect_aromatic_rings(
    atoms: List[Dict],
    bonds: List[Dict],
    atom_map: Dict[str, Dict],
    bonds_by_atom: Dict[str, List[Dict]]
) -> List[str]:
    """Detect aromatic rings (simplified - finds 6-membered rings)"""
    # Simple cycle detection
    visited = set()
    rings = []
    
    def find_cycle(start: str, path: List[str], target: str, depth: int):
        if depth > 6:
            return None
        if start == target and len(path) == 6:
            return path
        visited.add(start)
        for bond in bonds_by_atom.get(start, []):
            neighbor = bond["atom1"] if bond["atom1"] != start else bond["atom2"]
            if neighbor not in visited or (neighbor == target and len(path) == 5):
                result = find_cycle(neighbor, path + [neighbor], target, depth + 1)
                if result:
                    return result
        visited.remove(start)
        return None
    
    for atom in atoms:
        if atom["element"] == "C":
            cycle = find_cycle(atom["id"], [atom["id"]], atom["id"], 0)
            if cycle:
                rings.extend(cycle)
                break
    
    return list(set(rings)) if rings else []

