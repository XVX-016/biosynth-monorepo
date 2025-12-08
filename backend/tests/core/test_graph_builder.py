"""
Unit tests for MoleculeGraphBuilder.
"""
import pytest
from chem.core.graph_builder import MoleculeGraphBuilder


def test_from_molforge():
    """Test deserialization from .molforge format."""
    data = {
        "atoms": {
            "0": {"id": 0, "element": "C", "charge": 0},
            "1": {"id": 1, "element": "H", "charge": 0}
        },
        "positions": {
            "0": [0.0, 0.0, 0.0],
            "1": [1.0, 0.0, 0.0]
        },
        "bonds": [
            {"atom_a": 0, "atom_b": 1, "order": 1}
        ]
    }

    atoms, positions, bonds, adjacency = MoleculeGraphBuilder.from_molforge(data)

    assert len(atoms) == 2
    assert len(bonds) == 1
    assert atoms[0].element == "C"
    assert atoms[1].element == "H"
    assert positions[0] == (0.0, 0.0, 0.0)
    assert positions[1] == (1.0, 0.0, 0.0)
    assert bonds[0].atom_a == 0
    assert bonds[0].atom_b == 1
    assert bonds[0].order == 1
    assert 1 in adjacency[0]
    assert 0 in adjacency[1]
