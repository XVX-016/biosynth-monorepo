"""
Traditional chemistry rule validation engine.
"""
from typing import List, Dict
from ..core.models import Atom, Bond
from ..core.periodic_table import ELEMENTS


class RuleEngine:
    """Traditional chemistry constraints."""

    def validate_bond(self, atom_a: Atom, atom_b: Atom, order: int) -> bool:
        """
        Validate if a bond is chemically valid.
        
        Args:
            atom_a: First atom
            atom_b: Second atom
            order: Bond order (1=single, 2=double, 3=triple)
            
        Returns:
            True if bond is valid, False otherwise
        """
        # 1. Can't bond an atom to itself
        if atom_a.id == atom_b.id:
            return False

        # 2. Check if elements exist in periodic table
        if atom_a.element not in ELEMENTS or atom_b.element not in ELEMENTS:
            return False

        # 3. Maximum bond order rule (simplified)
        max_order = min(
            ELEMENTS[atom_a.element]["max_valence"],
            ELEMENTS[atom_b.element]["max_valence"]
        )
        if order > max_order:
            return False

        return True

    def filter_invalid(self, atoms: Dict[int, Atom], ml_bonds: List[Bond]) -> List[Bond]:
        """
        Drop ML bonds that violate chemistry rules.
        
        Args:
            atoms: Dictionary of atoms by ID
            ml_bonds: List of bonds predicted by ML
            
        Returns:
            List of valid bonds
        """
        valid = []
        for bond in ml_bonds:
            a = atoms.get(bond.atom_a)
            b = atoms.get(bond.atom_b)
            
            if a is None or b is None:
                continue
                
            if self.validate_bond(a, b, bond.order):
                valid.append(bond)

        return valid
