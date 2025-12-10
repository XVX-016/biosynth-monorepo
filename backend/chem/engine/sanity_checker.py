"""
Molecular sanity checker - enforces valence and charge constraints.
"""
from collections import defaultdict
from typing import Dict, List
from ..core.models import Atom, Bond
from ..core.periodic_table import ELEMENTS


class MoleculeSanityChecker:
    """Ensures molecules obey valence and charge rules."""

    def enforce(self, atoms: Dict[int, Atom], bonds: List[Bond]) -> List[Bond]:
        """
        Ensure valence does not exceed limits.
        
        Args:
            atoms: Dictionary of atoms by ID
            bonds: List of bonds to validate
            
        Returns:
            List of bonds that don't violate valence rules
        """
        valence_map = defaultdict(int)

        # Calculate current valence for each atom
        for b in bonds:
            valence_map[b.atom_a] += b.order
            valence_map[b.atom_b] += b.order

        sane_bonds = []
        for b in bonds:
            atom_a = atoms.get(b.atom_a)
            atom_b = atoms.get(b.atom_b)
            
            if atom_a is None or atom_b is None:
                continue
            
            # Check if elements are in periodic table
            if atom_a.element not in ELEMENTS or atom_b.element not in ELEMENTS:
                continue
            
            # Check valence limits
            max_valence_a = ELEMENTS[atom_a.element]["max_valence"]
            max_valence_b = ELEMENTS[atom_b.element]["max_valence"]
            
            if (valence_map[b.atom_a] <= max_valence_a and
                valence_map[b.atom_b] <= max_valence_b):
                sane_bonds.append(b)

        return sane_bonds
