"""
Molecule graph builder - converts serialized data to internal structures.
"""
import json
from typing import Dict, List, Tuple
from .models import Atom, Bond


class MoleculeGraphBuilder:
    """
    Converts raw serialized molecule JSON into
    Atom/Bond objects + normalized adjacency map.
    """

    @staticmethod
    def from_molforge(data: dict) -> Tuple[Dict[int, Atom], Dict[int, Tuple[float, float, float]], List[Bond], Dict[int, List[int]]]:
        """
        Parse .molforge JSON format.
        
        Returns:
            (atoms, positions, bonds, adjacency)
        """
        atoms = {
            int(k): Atom(id=v["id"], element=v["element"], charge=v.get("charge", 0))
            for k, v in data["atoms"].items()
        }

        positions = {
            int(k): tuple(v)
            for k, v in data["positions"].items()
        }

        bonds = [
            Bond(
                atom_a=b["atom_a"],
                atom_b=b["atom_b"],
                order=b["order"]
            )
            for b in data["bonds"]
        ]

        adjacency = {i: [] for i in atoms}
        for b in bonds:
            adjacency[b.atom_a].append(b.atom_b)
            adjacency[b.atom_b].append(b.atom_a)

        return atoms, positions, bonds, adjacency
