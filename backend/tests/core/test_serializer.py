"""
Unit tests for MoleculeSerializer.
"""
import json
import pytest
from chem.core.serializer import MoleculeSerializer
from chem.core.models import Atom, Bond


def test_molforge_roundtrip():
    """Test serialization to .molforge format."""
    atoms = {
        0: Atom(id=0, element="C"),
        1: Atom(id=1, element="H"),
    }
    positions = {0: (0.0, 0.0, 0.0), 1: (1.0, 0.0, 0.0)}
    bonds = [Bond(atom_a=0, atom_b=1, order=1)]

    data = MoleculeSerializer.to_molforge(atoms, bonds, positions)
    parsed = json.loads(data)

    assert "atoms" in parsed
    assert "bonds" in parsed
    assert "positions" in parsed
    assert len(parsed["atoms"]) == 2
    assert len(parsed["bonds"]) == 1
    assert parsed["atoms"]["0"]["element"] == "C"
    assert parsed["atoms"]["1"]["element"] == "H"


def test_mol_format():
    """Test serialization to .mol V2000 format."""
    atoms = {
        0: Atom(id=0, element="C"),
        1: Atom(id=1, element="H"),
    }
    positions = {
        0: (0.0, 0.0, 0.0),
        1: (1.0, 0.0, 0.0)
    }
    bonds = [Bond(atom_a=0, atom_b=1, order=1)]

    mol = MoleculeSerializer.to_mol(atoms, bonds, positions)
    assert "V2000" in mol
    assert "C" in mol
    assert "H" in mol
    assert "M  END" in mol


def test_pdb_format():
    """Test serialization to .pdb format."""
    atoms = {
        0: Atom(id=0, element="C"),
        1: Atom(id=1, element="H"),
    }
    positions = {
        0: (0.0, 0.0, 0.0),
        1: (1.0, 0.0, 0.0)
    }

    pdb = MoleculeSerializer.to_pdb(atoms, positions)
    assert "ATOM" in pdb
    assert "END" in pdb
    assert "C" in pdb
    assert "H" in pdb
