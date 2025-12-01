"""
Step 5: Generate Fingerprints for ML Models

Precomputes fingerprints using Phase 5 featurizers (ECFP, Morgan).
"""

import sys
from pathlib import Path
import csv
import json
import logging
import numpy as np

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

try:
    from backend.ai.featurizer import get_ecfp
    FEATURIZER_AVAILABLE = True
except ImportError:
    try:
        from src.ai.featurizer import get_ecfp
        FEATURIZER_AVAILABLE = True
    except ImportError:
        FEATURIZER_AVAILABLE = False
        logging.warning("Featurizer not available. Fingerprints will be mock.")

logger = logging.getLogger(__name__)


def generate_fingerprints():
    """
    Generate fingerprints for all molecules in the dataset.
    
    Uses Phase 5 featurizers to compute ECFP fingerprints.
    """
    data_dir = Path("data/datasets")
    cleaned_dir = data_dir / "cleaned"
    fingerprints_dir = data_dir / "fingerprints"
    fingerprints_dir.mkdir(parents=True, exist_ok=True)
    
    # Process train, validation, and test sets
    for split_name in ["train", "validation", "test"]:
        input_file = cleaned_dir / f"{split_name}.csv"
        
        if not input_file.exists():
            # Fallback to properties.csv
            input_file = data_dir / "properties.csv"
            if not input_file.exists():
                logger.warning(f"Input file not found for {split_name}")
                continue
        
        logger.info(f"Generating fingerprints for {split_name}...")
        
        fingerprints = []
        with open(input_file, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                smiles = row.get("smiles", "").strip()
                if not smiles:
                    continue
                
                # Generate fingerprint
                fp = _generate_fingerprint(smiles)
                if fp is not None:
                    fingerprints.append({
                        "smiles": smiles,
                        "fingerprint": fp,
                    })
        
        # Save fingerprints
        if fingerprints:
            output_file = fingerprints_dir / f"{split_name}_fingerprints.json"
            with open(output_file, "w") as f:
                json.dump(fingerprints, f, indent=2)
            
            logger.info(f"Generated {len(fingerprints)} fingerprints: {output_file}")
            
            # Also save as numpy array for faster loading
            np_file = fingerprints_dir / f"{split_name}_fingerprints.npy"
            fp_array = np.array([f["fingerprint"] for f in fingerprints])
            np.save(np_file, fp_array)
            
            # Save SMILES mapping
            smiles_file = fingerprints_dir / f"{split_name}_smiles.json"
            smiles_list = [f["smiles"] for f in fingerprints]
            with open(smiles_file, "w") as f:
                json.dump(smiles_list, f)
            
            logger.info(f"Saved numpy array: {np_file}")


def _generate_fingerprint(smiles: str) -> list:
    """
    Generate fingerprint for a SMILES string.
    
    Args:
        smiles: SMILES string
    
    Returns:
        Fingerprint as list of bit indices or None if failed
    """
    if not FEATURIZER_AVAILABLE:
        # Mock fingerprint
        return list(range(10))  # Mock: 10 bits
    
    try:
        fp_set = get_ecfp(smiles, radius=2, n_bits=2048)
        if fp_set:
            # Convert set to sorted list
            return sorted(list(fp_set))
        return None
    except Exception as e:
        logger.warning(f"Fingerprint generation error for {smiles}: {e}")
        return None


if __name__ == "__main__":
    generate_fingerprints()

