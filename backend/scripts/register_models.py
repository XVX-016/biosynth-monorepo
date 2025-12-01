"""
Step 5: Register Trained Models

Saves and registers trained models in Phase 5 model registry.
"""

import sys
from pathlib import Path
import logging
from typing import Dict, Any
from datetime import datetime

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

try:
    from ml.registry import ModelRegistry
    REGISTRY_AVAILABLE = True
except ImportError:
    REGISTRY_AVAILABLE = False
    logging.warning("Model registry not available.")

logger = logging.getLogger(__name__)


def register_models(
    trained_models: Dict[str, Any],
    evaluation_results: Dict[str, Dict[str, float]],
):
    """
    Register trained models in model registry.
    
    Args:
        trained_models: Dictionary of trained models
        evaluation_results: Evaluation metrics
    """
    if not REGISTRY_AVAILABLE:
        logger.warning("Model registry not available. Skipping registration.")
        return
    
    registry = ModelRegistry(registry_path="data/models/registry.json")
    
    for prop, model_info in trained_models.items():
        logger.info(f"Registering model for {prop}...")
        
        config = model_info["config"]
        metrics = evaluation_results.get(prop, {})
        
        # Create model ID
        model_id = f"attention-gnn-{prop.lower()}"
        
        # Register model
        registry.register(
            model_id=model_id,
            model_type="attention-gnn",
            path=model_info["path"],
            metadata={
                "node_feat_dim": config["node_feat_dim"],
                "edge_feat_dim": config["edge_feat_dim"],
                "hidden_dim": config["hidden_dim"],
                "out_dim": config["out_dim"],
                "num_layers": config["num_layers"],
                "heads": config["heads"],
                "property": prop,
                "training_date": datetime.now().isoformat(),
                "dataset": "prepared_datasets",
                "metrics": metrics,
                "best_val_loss": model_info.get("best_val_loss", 0.0),
            },
            description=f"Attention-GNN model for {prop} prediction",
            is_default=(prop == "logP"),  # Make logP default
        )
        
        logger.info(f"Registered {model_id}")
    
    logger.info("Model registration complete!")


if __name__ == "__main__":
    # This would be called from train_ml_models.py
    pass

