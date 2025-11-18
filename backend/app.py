"""
MolForge Backend API
FastAPI application entrypoint
"""
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from backend.config import settings
from backend.routes import predict as predict_router
from backend.routes import generate as generate_router
from backend.routes import library as library_router
from backend.routes import admin as admin_router
from backend.db import init_db
from backend.services.prediction_service import PredictionService
from backend.models.schemas.prediction_schema import PredictOut, PredictIn

app = FastAPI(
    title=settings.API_TITLE,
    version=settings.API_VERSION,
    description=settings.API_DESCRIPTION
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.CORS_ORIGINS,
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
    return {
        "message": "MolForge Backend API",
        "version": settings.API_VERSION
    }


@app.post("/predict-fast")
def predict_fast(payload: PredictFastIn):
    """
    Fast prediction using ONNX model
    """
    try:
        properties = PredictionService.predict_properties(payload.smiles, use_onnx=True)
        return PredictOut(properties=properties, smiles=payload.smiles)
    except Exception as e:
        # Fallback to regular predictor if ONNX fails
        try:
            properties = PredictionService.predict_properties(payload.smiles, use_onnx=False)
            return PredictOut(properties=properties, smiles=payload.smiles)
        except Exception as e2:
            raise HTTPException(status_code=500, detail=str(e2))


@app.get("/health")
def health():
    return {"status": "healthy"}

