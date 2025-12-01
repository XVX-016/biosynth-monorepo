"""
FastAPI router for molecular property prediction with attention

Integrates with PredictionEngine and ModelRegistry.
"""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List, Optional, Dict, Any
import logging

from ml.prediction_engine import PredictionEngine
from ml.registry import ModelRegistry

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api", tags=["prediction"])

# Initialize engines (singleton pattern)
_model_registry = None
_prediction_engine = None

def get_model_registry():
    global _model_registry
    if _model_registry is None:
        _model_registry = ModelRegistry(registry_path="data/models/registry.json")
    return _model_registry

def get_prediction_engine():
    global _prediction_engine
    if _prediction_engine is None:
        _prediction_engine = PredictionEngine(get_model_registry())
    return _prediction_engine


# Request/Response models
class PredictRequest(BaseModel):
    class Config:
        protected_namespaces = ()  # Fix Pydantic warning
    
    smiles: Optional[str] = None
    graph: Optional[Dict[str, Any]] = None
    molecule: Optional[Dict[str, Any]] = None
    node_order: Optional[List[str]] = None
    model_id: Optional[str] = None
    properties: Optional[List[str]] = None
    return_attention: bool = False
    attention_layer: str = "last"


class BatchPredictRequest(BaseModel):
    class Config:
        protected_namespaces = ()  # Fix Pydantic warning
    
    molecules: Optional[List[Dict[str, Any]]] = None  # List of molecule dicts with smiles or graph
    inputs: Optional[List[Dict[str, Any]]] = None  # Alternative format (for compatibility)
    model_id: Optional[str] = None
    properties: Optional[List[str]] = None
    batch_size: int = 32
    
    def __init__(self, **data):
        """Custom init to handle both molecules and inputs formats."""
        # If inputs is provided but molecules is not, copy inputs to molecules
        if 'inputs' in data and data.get('inputs') and 'molecules' not in data:
            data['molecules'] = data['inputs']
        super().__init__(**data)


class SimilarMoleculeRequest(BaseModel):
    smiles: str
    k: int = 5
    fingerprint_type: str = "ECFP4"


@router.post("/predict/property")
async def predict_property(request: PredictRequest):
    """
    Predict molecular properties.
    
    Accepts SMILES, graph JSON, or molecule JSON.
    Returns predictions with optional attention weights.
    """
    try:
        # Build input dict
        input_data = {}
        if request.smiles:
            input_data["smiles"] = request.smiles
            if request.node_order:
                input_data["node_order"] = request.node_order
        elif request.graph:
            input_data["graph"] = request.graph
        elif request.molecule:
            input_data["molecule"] = request.molecule
        else:
            raise HTTPException(status_code=400, detail="No molecule input provided")
        
        # Run prediction
        engine = get_prediction_engine()
        result = engine.predict(
            input_data=input_data,
            model_id=request.model_id,
            properties=request.properties,
            return_attention=request.return_attention,
            attention_layer=request.attention_layer,
        )
        
        return result.to_dict()
    
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Prediction error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Prediction failed: {str(e)}")


@router.post("/predict/all")
async def predict_all(request: PredictRequest):
    """
    Predict all available properties.
    """
    # Same as predict_property but with all properties
    request.properties = None  # Will use defaults
    return await predict_property(request)


@router.post("/predict/attention-map")
async def get_attention_map(request: PredictRequest):
    """
    Get attention map for molecule.
    
    Returns attention weights mapped to atoms/bonds.
    """
    request.return_attention = True
    
    result = await predict_property(request)
    
    # Format attention map
    attention_map = {
        "atoms": {},
        "bonds": {},
    }
    
    if "attentions" in result and "node_mapping" in result:
        node_mapping = result["node_mapping"]
        node_importance = result.get("attentions", {}).get("node_importance", [])
        
        # Map node importance to atom IDs
        for idx, importance in enumerate(node_importance):
            atom_id = node_mapping.get(idx)
            if atom_id:
                attention_map["atoms"][atom_id] = importance
        
        # Map edge attentions to bonds
        if "edge_index" in result and "attentions" in result:
            edge_index = result["edge_index"]
            edge_attentions = result["attentions"].get("layer_last", [])
            
            for e, attention in enumerate(edge_attentions):
                if e < len(edge_index[0]):
                    u_idx = edge_index[0][e]
                    v_idx = edge_index[1][e]
                    u_id = node_mapping.get(u_idx)
                    v_id = node_mapping.get(v_idx)
                    
                    if u_id and v_id:
                        # Create bond key (sorted for consistency)
                        bond_key = f"{min(u_id, v_id)}_{max(u_id, v_id)}"
                        attention_map["bonds"][bond_key] = attention
    
    return {
        "attention_map": attention_map,
        "predictions": result.get("predictions", {}),
        "model_id": result.get("model_id"),
    }


@router.post("/predict/batch")
async def predict_batch(request: BatchPredictRequest):
    """
    Batch prediction for multiple molecules.
    
    Accepts either:
    - molecules: List of molecule dicts with smiles or graph
    - inputs: Alternative format (for compatibility)
    """
    try:
        engine = get_prediction_engine()
        
        # Support both formats - __init__ already handles copying inputs to molecules
        inputs = request.molecules if request.molecules else []
        
        if not inputs:
            raise HTTPException(status_code=400, detail="No molecules provided. Use 'molecules' or 'inputs' field.")
        
        results = engine.predict_batch(
            inputs,
            model_id=request.model_id,
            batch_size=request.batch_size,
        )
        
        # Filter to requested properties if specified
        if request.properties:
            for result in results:
                # Only keep requested properties
                filtered_predictions = {k: v for k, v in result.predictions.items() if k in request.properties}
                result.predictions = filtered_predictions
        
        return {
            "predictions": [r.to_dict() for r in results],
            "count": len(results),
        }
    
    except Exception as e:
        logger.error(f"Batch prediction error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/suggest/similar")
async def suggest_similar(request: SimilarMoleculeRequest):
    """
    Find similar molecules using fingerprint similarity.
    
    TODO: Implement fingerprint search against database
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import DataStructs
        from rdkit.Chem.Fingerprints import FingerprintMols
        
        mol = Chem.MolFromSmiles(request.smiles)
        if mol is None:
            raise HTTPException(status_code=400, detail="Invalid SMILES")
        
        # Generate fingerprint
        fp = FingerprintMols.FingerprintMol(mol)
        
        # TODO: Search against database
        # For now, return placeholder
        return {
            "similar": [],
            "count": 0,
            "note": "Similarity search not yet implemented",
        }
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/ml/models")
async def list_models():
    """List all registered models."""
    registry = get_model_registry()
    models = registry.list_models()
    return {
        "models": models,
        "count": len(models),
    }


@router.get("/ml/models/active")
async def get_active_model():
    """Get currently active/default model."""
    registry = get_model_registry()
    engine = get_prediction_engine()
    model_id = registry.get_default_model_id()
    if model_id:
        model_info = registry.get_model(model_id)
        return {
            "model": model_info,
            "loaded": model_id in engine.loaded_models,
        }
    return {"model": None}

