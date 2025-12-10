"""
Unit tests for BondOptimizer.
"""
import pytest
from chem.engine.optimizer import BondOptimizer
from chem.core.models import Atom, Bond


def test_optimizer_passthrough():
    """Test that optimizer currently returns bonds unchanged."""
    opt = BondOptimizer()
    atoms = {0: Atom(id=0, element="C"), 1: Atom(id=1, element="H")}
    bonds = [Bond(atom_a=0, atom_b=1, order=1)]
    out = opt.optimize(atoms, bonds)
    assert out == bonds
    assert len(out) == 1
