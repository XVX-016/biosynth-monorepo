"""
Unit tests for BondEngine.
"""
import pytest
from chem.engine.bond_engine import BondEngine
from chem.core.models import Atom


def test_predict_basic_ch4():
    """Test bond prediction for methane (CH4)."""
    engine = BondEngine()

    atoms = {
        0: Atom(id=0, element="C"),
        1: Atom(id=1, element="H"),
        2: Atom(id=2, element="H"),
        3: Atom(id=3, element="H"),
        4: Atom(id=4, element="H"),
    }

    positions = {
        0: (0.0, 0.0, 0.0),
        1: (1.09, 0.0, 0.0),
        2: (-1.09, 0.0, 0.0),
        3: (0.0, 1.09, 0.0),
        4: (0.0, -1.09, 0.0),
    }

    bonds = engine.predict_all_bonds(atoms, positions)

    assert len(bonds) == 4
    for b in bonds:
        assert b.order == 1
        assert b.atom_a == 0 or b.atom_b == 0  # All bonds connect to carbon


def test_add_bond_valid():
    """Test manual bond creation with valid atoms."""
    engine = BondEngine()
    a = Atom(id=1, element="C")
    b = Atom(id=2, element="H")
    
    bond = engine.add_bond(a, b, 1)
    
    assert bond.atom_a == 1
    assert bond.atom_b == 2
    assert bond.order == 1


def test_add_bond_invalid():
    """Test that invalid bonds raise ValueError."""
    engine = BondEngine()
    a = Atom(id=1, element="C")
    
    # Can't bond atom to itself
    with pytest.raises(ValueError):
        engine.add_bond(a, a, 1)
