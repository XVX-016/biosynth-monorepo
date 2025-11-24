"""
Force field energy calculations (MMFF/GAFF-like)
"""
from typing import Dict, List, Any, Tuple
import math

# Force field parameters (simplified)
BOND_PARAMS = {
    ("C", "C"): {"k": 350.0, "r0": 1.54},
    ("C", "H"): {"k": 340.0, "r0": 1.09},
    ("C", "O"): {"k": 320.0, "r0": 1.43},
    ("C", "N"): {"k": 310.0, "r0": 1.47},
    ("O", "H"): {"k": 450.0, "r0": 0.96},
    ("N", "H"): {"k": 430.0, "r0": 1.01},
}

ANGLE_PARAMS = {
    ("C", "C", "C"): {"k": 40.0, "theta0": 109.5},
    ("C", "C", "H"): {"k": 35.0, "theta0": 109.5},
    ("C", "O", "H"): {"k": 50.0, "theta0": 104.5},
}

DIHEDRAL_PARAMS = {
    ("C", "C", "C", "C"): {"v": [0.0, 0.0, 0.5, 0.0]},  # Simplified
}

# Lennard-Jones parameters (ε in kcal/mol, σ in Å)
LJ_PARAMS = {
    "C": {"epsilon": 0.1094, "sigma": 3.399},
    "H": {"epsilon": 0.0157, "sigma": 2.571},
    "O": {"epsilon": 0.2100, "sigma": 3.066},
    "N": {"epsilon": 0.1700, "sigma": 3.250},
}

def distance(atom1: Dict, atom2: Dict) -> float:
    """Calculate distance between two atoms"""
    pos1 = atom1["position"]
    pos2 = atom2["position"]
    dx = pos2[0] - pos1[0]
    dy = pos2[1] - pos1[1]
    dz = pos2[2] - pos1[2]
    return math.sqrt(dx*dx + dy*dy + dz*dz)

def angle(atom1: Dict, atom2: Dict, atom3: Dict) -> float:
    """Calculate angle between three atoms (in degrees)"""
    v1 = [
        atom1["position"][0] - atom2["position"][0],
        atom1["position"][1] - atom2["position"][1],
        atom1["position"][2] - atom2["position"][2]
    ]
    v2 = [
        atom3["position"][0] - atom2["position"][0],
        atom3["position"][1] - atom2["position"][1],
        atom3["position"][2] - atom2["position"][2]
    ]
    
    dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
    len1 = math.sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2])
    len2 = math.sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2])
    
    cos_theta = dot / (len1 * len2) if len1 * len2 > 0 else 0
    cos_theta = max(-1, min(1, cos_theta))  # Clamp
    return math.degrees(math.acos(cos_theta))

def calculate_bond_energy(molecule: Dict[str, Any]) -> Tuple[float, Dict[str, float]]:
    """Calculate bond stretching energy"""
    atoms = {a["id"]: a for a in molecule["atoms"]}
    bonds = molecule["bonds"]
    
    total_energy = 0.0
    per_bond = {}
    
    for bond in bonds:
        a1 = atoms.get(bond["atom1"])
        a2 = atoms.get(bond["atom2"])
        if not a1 or not a2:
            continue
        
        elem1 = a1["element"]
        elem2 = a2["element"]
        key = tuple(sorted([elem1, elem2]))
        
        params = BOND_PARAMS.get(key) or BOND_PARAMS.get(("C", "C"))
        if not params:
            continue
        
        r = distance(a1, a2)
        r0 = params["r0"]
        k = params["k"]
        
        # Harmonic potential: E = 0.5 * k * (r - r0)^2
        energy = 0.5 * k * (r - r0) ** 2
        total_energy += energy
        per_bond[bond["id"]] = energy
    
    return total_energy, per_bond

def calculate_angle_energy(molecule: Dict[str, Any]) -> Tuple[float, Dict[str, float]]:
    """Calculate angle bending energy"""
    atoms = {a["id"]: a for a in molecule["atoms"]}
    bonds = molecule["bonds"]
    
    # Build adjacency
    adj = {a["id"]: [] for a in molecule["atoms"]}
    for bond in bonds:
        adj[bond["atom1"]].append(bond["atom2"])
        adj[bond["atom2"]].append(bond["atom1"])
    
    total_energy = 0.0
    per_angle = {}
    angle_id = 0
    
    for atom_id, neighbors in adj.items():
        if len(neighbors) < 2:
            continue
        
        atom = atoms[atom_id]
        for i in range(len(neighbors)):
            for j in range(i + 1, len(neighbors)):
                a1 = atoms[neighbors[i]]
                a2 = atom
                a3 = atoms[neighbors[j]]
                
                key = tuple(sorted([a1["element"], a2["element"], a3["element"]]))
                params = ANGLE_PARAMS.get(key) or ANGLE_PARAMS.get(("C", "C", "C"))
                if not params:
                    continue
                
                theta = angle(a1, a2, a3)
                theta0 = params["theta0"]
                k = params["k"]
                
                # Harmonic potential: E = 0.5 * k * (θ - θ0)^2
                dtheta = math.radians(theta - theta0)
                energy = 0.5 * k * dtheta ** 2
                total_energy += energy
                per_angle[f"angle_{angle_id}"] = energy
                angle_id += 1
    
    return total_energy, per_angle

def calculate_nonbonded_energy(molecule: Dict[str, Any]) -> float:
    """Calculate Lennard-Jones + Coulomb nonbonded energy"""
    atoms = molecule["atoms"]
    bonds = molecule["bonds"]
    
    # Build bonded pairs set
    bonded_pairs = set()
    for bond in bonds:
        bonded_pairs.add(tuple(sorted([bond["atom1"], bond["atom2"]])))
    
    total_energy = 0.0
    
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            a1 = atoms[i]
            a2 = atoms[j]
            
            # Skip bonded pairs
            pair_key = tuple(sorted([a1["id"], a2["id"]]))
            if pair_key in bonded_pairs:
                continue
            
            r = distance(a1, a2)
            if r < 0.1:  # Avoid division by zero
                continue
            
            # Lennard-Jones
            elem1 = a1["element"]
            elem2 = a2["element"]
            lj1 = LJ_PARAMS.get(elem1, LJ_PARAMS["C"])
            lj2 = LJ_PARAMS.get(elem2, LJ_PARAMS["C"])
            
            epsilon = math.sqrt(lj1["epsilon"] * lj2["epsilon"])
            sigma = (lj1["sigma"] + lj2["sigma"]) / 2
            
            sr = sigma / r
            sr6 = sr ** 6
            sr12 = sr6 ** 2
            
            lj_energy = 4 * epsilon * (sr12 - sr6)
            total_energy += lj_energy
    
    return total_energy

def calculate_total_energy(molecule: Dict[str, Any]) -> Dict[str, Any]:
    """
    Calculate total molecular mechanics energy
    
    Returns:
        {
            "total_energy": float,
            "bond_energy": float,
            "angle_energy": float,
            "nonbonded_energy": float,
            "per_bond": {...},
            "per_angle": {...}
        }
    """
    bond_energy, per_bond = calculate_bond_energy(molecule)
    angle_energy, per_angle = calculate_angle_energy(molecule)
    nonbonded_energy = calculate_nonbonded_energy(molecule)
    
    total_energy = bond_energy + angle_energy + nonbonded_energy
    
    return {
        "total_energy": total_energy,
        "bond_energy": bond_energy,
        "angle_energy": angle_energy,
        "nonbonded_energy": nonbonded_energy,
        "per_bond": per_bond,
        "per_angle": per_angle
    }

