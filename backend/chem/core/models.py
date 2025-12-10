"""
Core molecular data models for MolForge.
"""
from dataclasses import dataclass
from typing import Optional


@dataclass
class Atom:
    """Represents an atom in a molecule."""
    id: int
    element: str
    charge: int = 0
    x: float = 0.0
    y: float = 0.0
    z: float = 0.0

    def to_dict(self):
        return {
            "id": self.id,
            "element": self.element,
            "charge": self.charge,
            "x": self.x,
            "y": self.y,
            "z": self.z
        }


@dataclass
class Bond:
    """Represents a bond between two atoms."""
    atom_a: int
    atom_b: int
    order: int = 1

    def to_dict(self):
        return {
            "atom_a": self.atom_a,
            "atom_b": self.atom_b,
            "order": self.order
        }
