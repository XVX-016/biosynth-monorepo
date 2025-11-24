"""
MolForge Backend API
FastAPI application entrypoint
"""
from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from pydantic import BaseModel
import traceback
from backend.config import settings
from backend.routes import predict as predict_router
from backend.routes import generate as generate_router
from backend.routes import library as library_router
from backend.routes import admin as admin_router
from backend.routes import convert as convert_router
from backend.routes import thumbnails as thumbnails_router
from backend.routes import relax as relax_router
from backend.db import init_db
from backend.services.prediction_service import PredictionService
from backend.models.schemas.prediction_schema import PredictOut, PredictIn

app = FastAPI(
    title=settings.API_TITLE,
    version=settings.API_VERSION,
    description=settings.API_DESCRIPTION
)

# CORS middleware must be added before routes
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.CORS_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
    expose_headers=["*"],
)

# Global exception handler to ensure CORS headers on errors
@app.exception_handler(Exception)
async def global_exception_handler(request: Request, exc: Exception):
    """Handle all exceptions and ensure CORS headers are included"""
    import traceback
    error_detail = str(exc)
    if settings.LOG_LEVEL == "DEBUG":
        error_detail += f"\n{traceback.format_exc()}"
    
    return JSONResponse(
        status_code=500,
        content={
            "detail": error_detail,
            "path": str(request.url),
        },
        headers={
            "Access-Control-Allow-Origin": request.headers.get("origin", "*"),
            "Access-Control-Allow-Credentials": "true",
        }
    )

# Initialize database
init_db()

# Mount routers
app.include_router(predict_router.router, prefix="/predict", tags=["predict"])
app.include_router(generate_router.router, prefix="/generate", tags=["generate"])
app.include_router(library_router.router, tags=["molecules"])
app.include_router(admin_router.router, tags=["admin"])
app.include_router(convert_router.router, tags=["convert"])
app.include_router(thumbnails_router.router, tags=["thumbnails"])
app.include_router(relax_router.router, prefix="/api", tags=["relax"])


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

