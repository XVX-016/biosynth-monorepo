"""
Molecule tokenizer - converts molecule to graph representation and fingerprints
"""
from typing import Dict, List, Set, Any
import hashlib

def tokenize_molecule(molecule: Dict[str, Any]) -> Dict[str, Any]:
    """
    Convert molecule (atoms + bonds) into graph representation and fingerprint
    
    Args:
        molecule: { "atoms": [...], "bonds": [...] }
    
    Returns:
        {
            "atoms": [...],
            "bonds": [...],
            "graph": {node_id: [neighbor_ids]},
            "fingerprint": bit_vector_string
        }
    """
    atoms = molecule.get("atoms", [])
    bonds = molecule.get("bonds", [])
    
    # Build adjacency graph
    graph: Dict[str, List[str]] = {}
    for atom in atoms:
        graph[atom["id"]] = []
    
    for bond in bonds:
        a1 = bond["atom1"]
        a2 = bond["atom2"]
        if a1 in graph:
            graph[a1].append(a2)
        if a2 in graph:
            graph[a2].append(a1)
    
    # Generate simple fingerprint (path-based)
    fingerprint = generate_fingerprint(atoms, bonds)
    
    return {
        "atoms": atoms,
        "bonds": bonds,
        "graph": graph,
        "fingerprint": fingerprint
    }

def generate_fingerprint(atoms: List[Dict], bonds: List[Dict]) -> str:
    """
    Generate a simple path-based fingerprint
    Returns a hash string representing the molecule structure
    """
    # Create canonical representation
    atom_strs = sorted([f"{a['element']}:{a['id']}" for a in atoms])
    bond_strs = sorted([
        f"{b['atom1']}-{b['atom2']}:{b.get('order', 1)}" 
        for b in bonds
    ])
    
    canonical = "|".join(atom_strs) + "||" + "|".join(bond_strs)
    return hashlib.md5(canonical.encode()).hexdigest()

