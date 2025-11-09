# backend/test_library.py
from fastapi.testclient import TestClient
from app import app

client = TestClient(app)

def test_save_and_get_molecule():
    payload = {
        "name": "TestMol",
        "smiles": "CCO",
        "json_graph": '{"atoms":[{"id":"atom_1","element":"C","position":[0,0,0]}],"bonds":[]}',
        "coords": "{}",
        "properties": "{}",
        "thumbnail_b64": None
    }
    r = client.post("/molecules/save", json=payload)
    assert r.status_code == 200
    id = r.json().get("id")
    assert id is not None

    r2 = client.get(f"/molecules/{id}")
    assert r2.status_code == 200
    assert r2.json()["name"] == "TestMol"

    r3 = client.delete(f"/molecules/{id}")
    assert r3.status_code == 200

def test_list_molecules():
    r = client.get("/molecules/list")
    assert r.status_code == 200
    assert isinstance(r.json(), list)

