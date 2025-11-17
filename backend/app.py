from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from backend.routes import predict as predict_router
from backend.routes import generate as generate_router
from backend.routes import library as library_router
from backend.routes import admin as admin_router
from backend.db import init_db

app = FastAPI(title="BioSynth AI Backend", version="0.1.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://localhost:5174", "http://127.0.0.1:5173", "http://127.0.0.1:5174"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Initialize database
init_db()

# Mount routers
app.include_router(predict_router.router, prefix="/predict", tags=["predict"])
app.include_router(generate_router.router, prefix="/generate", tags=["generate"])
app.include_router(library_router.router, tags=["molecules"])
app.include_router(admin_router.router, tags=["admin"])


class PredictFastIn(BaseModel):
    smiles: str


@app.get("/")
def root():
    return {"message": "BioSynth AI Backend API", "version": "0.1.0"}


@app.post("/predict-fast")
def predict_fast(payload: PredictFastIn):
    """
    Fast prediction using ONNX model
    """
    from backend.models.onnx_predictor import get_onnx_predictor
    from backend.utils.featurizer import featurize_smiles, validate_smiles
    from backend.routes.predict import PredictOut
    
    if not validate_smiles(payload.smiles):
        raise HTTPException(status_code=400, detail="Invalid SMILES string")
    
    features = featurize_smiles(payload.smiles)
    if features is None:
        raise HTTPException(status_code=400, detail="Failed to featurize SMILES")
    
    try:
        predictor = get_onnx_predictor()
        properties = predictor.predict(features)
        return PredictOut(properties=properties)
    except FileNotFoundError as e:
        # Fallback to regular predictor if ONNX not available
        from backend.routes.predict import predict
        from backend.routes.predict import PredictIn
        return predict(PredictIn(smiles=payload.smiles))


@app.get("/health")
def health():
    return {"status": "healthy"}

