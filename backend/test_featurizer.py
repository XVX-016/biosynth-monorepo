import pytest
from backend.ai.featurizer import featurize_smiles, validate_smiles


def test_featurize_valid_smiles():
    """Test featurization of valid SMILES"""
    features = featurize_smiles("CCO")
    assert features is not None
    assert len(features) == 2048
    assert all(isinstance(x, float) for x in features)


def test_featurize_invalid_smiles():
    """Test featurization of invalid SMILES"""
    features = featurize_smiles("INVALID")
    assert features is None


def test_validate_smiles_valid():
    """Test validation of valid SMILES"""
    assert validate_smiles("CCO") is True
    assert validate_smiles("C") is True


def test_validate_smiles_invalid():
    """Test validation of invalid SMILES"""
    assert validate_smiles("INVALID") is False
    assert validate_smiles("") is False

