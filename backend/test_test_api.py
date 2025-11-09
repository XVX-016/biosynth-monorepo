import pytest
from fastapi.testclient import TestClient
import app

client = TestClient(app.app)

def test_root():
    resp = client.get("/")
    assert resp.status_code == 200
    assert "message" in resp.json()

def test_health():
    resp = client.get("/health")
    assert resp.status_code == 200
    assert resp.json()["status"] == "healthy"

def test_predict_valid():
    resp = client.post("/predict", json={"smiles": "CCO"})
    assert resp.status_code == 200
    assert "properties" in resp.json()
    props = resp.json()["properties"]
    assert "stability" in props
    assert "toxicity" in props

def test_predict_invalid_smiles():
    resp = client.post("/predict", json={"smiles": "INVALID"})
    assert resp.status_code == 200
    assert "error" in resp.json()


