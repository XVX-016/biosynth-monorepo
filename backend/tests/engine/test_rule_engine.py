"""
Unit tests for RuleEngine.
"""
import pytest
from chem.engine.rule_engine import RuleEngine
from chem.core.models import Atom, Bond


def test_rule_engine_valid():
    """Test validation of valid C-H bond."""
    re = RuleEngine()
    a = Atom(id=1, element="C")
    b = Atom(id=2, element="H")
    assert re.validate_bond(a, b, 1) is True


def test_rule_engine_invalid_self_bond():
    """Test that self-bonds are rejected."""
    re = RuleEngine()
    a = Atom(id=1, element="C")
    assert re.validate_bond(a, a, 1) is False


def test_rule_engine_invalid_order():
    """Test that excessive bond orders are rejected."""
    re = RuleEngine()
    a = Atom(id=1, element="H")
    b = Atom(id=2, element="H")
    # H can only form single bonds
    assert re.validate_bond(a, b, 2) is False


def test_filter_invalid():
    """Test filtering of invalid bonds."""
    re = RuleEngine()
    atoms = {
        0: Atom(id=0, element="C"),
        1: Atom(id=1, element="H"),
    }
    
    bonds = [
        Bond(atom_a=0, atom_b=1, order=1),  # Valid
        Bond(atom_a=0, atom_b=0, order=1),  # Invalid (self-bond)
    ]
    
    valid = re.filter_invalid(atoms, bonds)
    assert len(valid) == 1
    assert valid[0].atom_a == 0
    assert valid[0].atom_b == 1
