import pytest
from fastapi.testclient import TestClient
from backend.app import app

client = TestClient(app)


def test_predict_valid():
    """Test prediction with valid SMILES"""
    resp = client.post("/predict/", json={"smiles": "CCO"})
    assert resp.status_code == 200
    assert "properties" in resp.json()
    props = resp.json()["properties"]
    assert "stability" in props
    assert "toxicity" in props
    assert "solubility" in props
    assert "bioavailability" in props
    assert "novelty" in props


def test_predict_invalid_smiles():
    """Test prediction with invalid SMILES"""
    resp = client.post("/predict/", json={"smiles": "INVALID"})
    assert resp.status_code == 400
    assert "detail" in resp.json()


def test_generate():
    """Test molecule generation"""
    resp = client.post("/generate/", json={"prompt": "Generate a molecule"})
    assert resp.status_code == 200
    assert "smiles" in resp.json()


def test_predict_fast():
    """Test fast prediction endpoint"""
    resp = client.post("/predict-fast", json={"smiles": "CCO"})
    assert resp.status_code == 200
    assert "properties" in resp.json()

