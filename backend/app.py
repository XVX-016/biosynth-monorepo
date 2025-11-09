from fastapi import FastAPI
from pydantic import BaseModel
from backend.routes import predict as predict_router
from backend.routes import generate as generate_router

app = FastAPI(title="BioSynth AI Backend", version="0.1.0")

# Mount routers
app.include_router(predict_router.router, prefix="/predict", tags=["predict"])
app.include_router(generate_router.router, prefix="/generate", tags=["generate"])


class PredictFastIn(BaseModel):
    smiles: str


@app.get("/")
def root():
    return {"message": "BioSynth AI Backend API", "version": "0.1.0"}


@app.post("/predict-fast")
def predict_fast(payload: PredictFastIn):
    """
    Fast prediction using ONNX model
    TODO: Implement ONNX inference
    """
    # Placeholder: use regular predict for now
    from backend.routes.predict import predict
    from backend.routes.predict import PredictIn
    
    return predict(PredictIn(smiles=payload.smiles))


@app.get("/health")
def health():
    return {"status": "healthy"}

