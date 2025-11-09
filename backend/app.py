from fastapi import FastAPI
from pydantic import BaseModel
import torch
from rdkit import Chem
from rdkit.Chem import AllChem
import os

app = FastAPI(title="BioSynth AI Backend", version="0.1.0")

class PredictIn(BaseModel):
    smiles: str

# Dummy model - replace with real PyTorch model loaded from weights
class DummyModel:
    def __call__(self, x):
        # Returns fake properties: [stability, toxicity, solubility, bioavailability, novelty]
        return [0.1, 0.2, 0.3, 0.4, 0.5]

model = DummyModel()

def featurize(smiles: str):
    """Convert SMILES string to Morgan fingerprint."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    arr = list(fp.ToBitString())
    return [float(x) for x in arr]

@app.get("/")
def root():
    return {"message": "BioSynth AI Backend API", "version": "0.1.0"}

@app.post("/predict")
def predict(payload: PredictIn):
    """Predict molecular properties from SMILES string."""
    feats = featurize(payload.smiles)
    if feats is None:
        return {"error": "invalid smiles"}
    # Convert into tensor (placeholder)
    # x = torch.tensor(feats).float().unsqueeze(0)
    # y = model(x)
    y = model(feats)
    return {
        "properties": {
            "stability": y[0],
            "toxicity": y[1],
            "solubility": y[2],
            "bioavailability": y[3],
            "novelty": y[4],
        }
    }

@app.get("/health")
def health():
    return {"status": "healthy"}

