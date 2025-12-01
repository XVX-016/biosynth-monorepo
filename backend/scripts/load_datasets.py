"""
Step 1: Load Prepared Datasets

Loads standardized training, validation, and test datasets
along with associated fingerprints.
"""

import sys
from pathlib import Path
import csv
import json
import numpy as np
import logging
from typing import Dict, List, Tuple, Any

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

logger = logging.getLogger(__name__)


def load_datasets() -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
    """
    Load training, validation, and test datasets.
    
    Returns:
        Tuple of (train_data, val_data, test_data) dictionaries
    """
    data_dir = Path("data/datasets")
    cleaned_dir = data_dir / "cleaned"
    fingerprints_dir = data_dir / "fingerprints"
    
    logger.info("Loading datasets from prepared files...")
    
    # Load training set
    train_data = _load_split("train", cleaned_dir, fingerprints_dir)
    logger.info(f"Training set: {len(train_data['smiles'])} molecules")
    
    # Load validation set
    val_data = _load_split("validation", cleaned_dir, fingerprints_dir)
    logger.info(f"Validation set: {len(val_data['smiles'])} molecules")
    
    # Load test set
    test_data = _load_split("test", cleaned_dir, fingerprints_dir)
    logger.info(f"Test set: {len(test_data['smiles'])} molecules")
    
    return train_data, val_data, test_data


def _load_split(
    split_name: str,
    cleaned_dir: Path,
    fingerprints_dir: Path,
) -> Dict[str, Any]:
    """
    Load a single split (train/validation/test).
    
    Args:
        split_name: "train", "validation", or "test"
        cleaned_dir: Directory with cleaned CSV files
        fingerprints_dir: Directory with fingerprint files
    
    Returns:
        Dictionary with smiles, properties, and fingerprints
    """
    csv_file = cleaned_dir / f"{split_name}.csv"
    
    if not csv_file.exists():
        raise FileNotFoundError(f"Dataset file not found: {csv_file}")
    
    smiles = []
    properties = {
        "logP": [],
        "solubility": [],
        "toxicity": [],
        "molecular_weight": [],
    }
    
    # Load CSV
    with open(csv_file, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            smi = row.get("smiles", "").strip()
            if not smi:
                continue
            
            smiles.append(smi)
            
            # Extract properties
            for prop in properties.keys():
                value = row.get(prop, "")
                try:
                    properties[prop].append(float(value) if value else 0.0)
                except (ValueError, TypeError):
                    properties[prop].append(0.0)
    
    # Load fingerprints if available
    fingerprints = None
    fp_file = fingerprints_dir / f"{split_name}_fingerprints.json"
    if fp_file.exists():
        try:
            with open(fp_file, "r") as f:
                fp_data = json.load(f)
            # Extract fingerprints (assuming format: [{"smiles": "...", "fingerprint": [...]}])
            fingerprints = []
            for entry in fp_data:
                if entry.get("smiles") in smiles:
                    fingerprints.append(entry.get("fingerprint", []))
            
            # Pad to match smiles length
            while len(fingerprints) < len(smiles):
                fingerprints.append([])
        except Exception as e:
            logger.warning(f"Failed to load fingerprints: {e}")
    
    return {
        "smiles": smiles,
        "properties": properties,
        "fingerprints": fingerprints,
        "count": len(smiles),
    }


if __name__ == "__main__":
    train, val, test = load_datasets()
    print(f"Train: {train['count']}, Val: {val['count']}, Test: {test['count']}")

