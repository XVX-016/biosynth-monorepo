"""
Unit tests for MoleculeSanityChecker.
"""
import pytest
from chem.engine.sanity_checker import MoleculeSanityChecker
from chem.core.models import Atom, Bond


def test_excess_valence_removed():
    """Test that bonds exceeding valence are removed."""
    checker = MoleculeSanityChecker()
    atoms = {
        0: Atom(id=0, element="C"),
        1: Atom(id=1, element="C")
    }

    bonds = [
        Bond(atom_a=0, atom_b=1, order=3),
        Bond(atom_a=0, atom_b=1, order=2),  # Illegal — total 5 > valence limit (4)
    ]

    sane = checker.enforce(atoms, bonds)
    # Both bonds should be removed since total valence (5) exceeds C limit (4)
    assert len(sane) == 0


def test_valid_valence_preserved():
    """Test that valid bonds are preserved."""
    checker = MoleculeSanityChecker()
    atoms = {
        0: Atom(id=0, element="C"),
        1: Atom(id=1, element="H"),
        2: Atom(id=2, element="H"),
    }

    bonds = [
        Bond(atom_a=0, atom_b=1, order=1),
        Bond(atom_a=0, atom_b=2, order=1),
    ]

    sane = checker.enforce(atoms, bonds)
    assert len(sane) == 2
