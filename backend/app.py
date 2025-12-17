"""
MolForge Backend API
FastAPI application entrypoint
"""
import sys
from pathlib import Path

# Ensure backend root is in Python path so 'chem' and other modules are importable
# This allows imports like 'from chem.search.substructure import find_substructure'
ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from pydantic import BaseModel
import traceback
from config import settings
from backend.routes import predict as predict_router
from backend.routes import generate as generate_router
from backend.routes import library as library_router
from backend.routes import admin as admin_router
from backend.routes import convert as convert_router
from backend.routes import thumbnails as thumbnails_router
from backend.routes import relax as relax_router
from backend.routes import search as search_router
from backend.routes import spectroscopy as spectroscopy_router
from backend.routes import energy as energy_router
from backend.routes import reaction as reaction_router
from backend.routes import retrosynthesis as retrosynthesis_router
from backend.routes import kab as kab_router
from backend.routes import quantum as quantum_router
from backend.routes import ml as ml_router
from backend.routes import collaboration as collaboration_router
from backend.routes import dashboard as dashboard_router
from backend.lab import routes as lab_router
from backend.api import predict as predict_api_router
from backend.routes import screening as screening_router
from backend.routes import search_phase7 as search_phase7_router
from backend.routes import qm_md as qm_md_router
from backend.api import search as search_api_router
from backend.api import screening as screening_api_router
from backend.api import conformers as conformers_api_router
from backend.api import orchestrator as orchestrator_api_router
from backend.api import phase10 as phase10_api_router
from backend.api import molecule as molecule_api_router
from backend.api import mentor as mentor_api_router
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
app.include_router(search_router.router, prefix="/api/search", tags=["search"])
app.include_router(spectroscopy_router.router, prefix="/api/spectroscopy", tags=["spectroscopy"])
app.include_router(energy_router.router, prefix="/api/energy", tags=["energy"])
app.include_router(reaction_router.router, prefix="/api/reaction", tags=["reaction"])
app.include_router(retrosynthesis_router.router, prefix="/api/retrosynthesis", tags=["retrosynthesis"])
app.include_router(kab_router.router, prefix="/api/kab", tags=["kab"])
app.include_router(quantum_router.router, prefix="/api/quantum", tags=["quantum"])
app.include_router(ml_router.router, prefix="/api/ml", tags=["ml"])
app.include_router(predict_api_router.router, tags=["prediction"])
app.include_router(screening_router.router, prefix="/api/screening", tags=["screening"])
app.include_router(search_phase7_router.router, prefix="/api/search", tags=["search-phase7"])
app.include_router(search_api_router.router, prefix="/api/search", tags=["search"])
app.include_router(screening_api_router.router, prefix="/api/screening", tags=["screening"])
app.include_router(conformers_api_router.router, prefix="/api/conformers", tags=["conformers"])
app.include_router(qm_md_router.router, prefix="/api", tags=["qm-md"])
app.include_router(orchestrator_api_router.router, prefix="/api/orchestrator", tags=["orchestrator"])
app.include_router(phase10_api_router.router, prefix="/api/phase10", tags=["phase10"])
app.include_router(molecule_api_router.router, tags=["molecule"])
app.include_router(collaboration_router.router, prefix="/api/collaboration", tags=["collaboration"])
app.include_router(dashboard_router.router, prefix="/api/dashboard", tags=["dashboard"])
app.include_router(lab_router.router, prefix="/api/lab", tags=["lab"])
app.include_router(mentor_api_router.router, prefix="/api/mentor", tags=["mentor"])


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

