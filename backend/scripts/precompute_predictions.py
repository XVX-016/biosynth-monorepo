"""
Step 7: Precompute Predictions for Screening

Precomputes property predictions for library molecules.
"""

import sys
from pathlib import Path
import csv
import logging
from typing import List, Dict, Any

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

try:
    from ml.prediction_engine import PredictionEngine
    from ml.registry import ModelRegistry
    PREDICTION_AVAILABLE = True
except ImportError:
    PREDICTION_AVAILABLE = False
    logging.warning("Prediction engine not available.")

logger = logging.getLogger(__name__)


def precompute_predictions():
    """
    Precompute property predictions for library molecules.
    
    Reads from /data/libraries/*.smi and stores results in
    /data/libraries/predictions.csv
    """
    if not PREDICTION_AVAILABLE:
        logger.warning("Prediction engine not available. Skipping precomputation.")
        return
    
    libraries_dir = Path("data/libraries")
    output_file = libraries_dir / "predictions.csv"
    
    # Find .smi files
    smi_files = list(libraries_dir.glob("*.smi"))
    if not smi_files:
        logger.warning(f"No .smi files found in {libraries_dir}")
        return
    
    logger.info(f"Precomputing predictions for library molecules...")
    
    # Initialize prediction engine
    registry = ModelRegistry()
    engine = PredictionEngine(registry)
    
    # Properties to predict
    properties = ["logP", "solubility", "toxicity"]
    
    # Read molecules from .smi files
    molecules = []
    for smi_file in smi_files:
        logger.info(f"Reading {smi_file}...")
        with open(smi_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split()
                if len(parts) >= 1:
                    smiles = parts[0]
                    name = parts[1] if len(parts) > 1 else smiles
                    molecules.append({"smiles": smiles, "name": name})
    
    logger.info(f"Found {len(molecules)} molecules")
    
    # Predict properties
    predictions = []
    for i, mol in enumerate(molecules):
        if (i + 1) % 10 == 0:
            logger.info(f"Processing {i+1}/{len(molecules)}...")
        
        try:
            result = engine.predict(
                input_data={"smiles": mol["smiles"]},
                properties=properties,
            )
            
            pred_dict = {
                "smiles": mol["smiles"],
                "name": mol["name"],
            }
            
            # Add predictions
            for prop in properties:
                pred_dict[prop] = result.predictions.get(prop, 0.0)
            
            predictions.append(pred_dict)
        except Exception as e:
            logger.warning(f"Failed to predict for {mol['smiles']}: {e}")
    
    # Save predictions
    if predictions:
        fieldnames = ["smiles", "name"] + properties
        with open(output_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(predictions)
        
        logger.info(f"Saved {len(predictions)} predictions to {output_file}")
    else:
        logger.error("No predictions generated!")


if __name__ == "__main__":
    precompute_predictions()

