"""
Step 1: Collect Property Datasets

Downloads and extracts SMILES strings with property values from
ChEMBL, ZINC, or other sources.
"""

import sys
from pathlib import Path
import csv
import logging

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available. SMILES validation will be limited.")

try:
    from backend.chem.utils.validators import validate_smiles
except ImportError:
    # Fallback: simple validation
    def validate_smiles(smiles: str) -> bool:
        return bool(smiles and len(smiles) > 0)

logger = logging.getLogger(__name__)


def collect_property_datasets():
    """
    Collect property datasets from various sources.
    
    For now, generates a sample dataset. In production, this would:
    1. Download from ChEMBL API
    2. Download from ZINC database
    3. Extract from custom sources
    """
    data_dir = Path("data/datasets")
    data_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = data_dir / "properties.csv"
    
    logger.info(f"Collecting property datasets to {output_file}")
    
    # Sample molecules with properties
    # In production, this would fetch from ChEMBL/ZINC APIs
    sample_data = [
        {
            "smiles": "CCO",
            "name": "ethanol",
            "logP": -0.31,
            "solubility": -0.77,  # -logS
            "toxicity": 0.0,
            "molecular_weight": 46.07,
        },
        {
            "smiles": "CC(=O)O",
            "name": "acetic_acid",
            "logP": -0.17,
            "solubility": 0.18,
            "toxicity": 0.0,
            "molecular_weight": 60.05,
        },
        {
            "smiles": "c1ccccc1",
            "name": "benzene",
            "logP": 2.13,
            "solubility": 1.64,
            "toxicity": 0.8,
            "molecular_weight": 78.11,
        },
        {
            "smiles": "CCc1ccccc1",
            "name": "toluene",
            "logP": 2.73,
            "solubility": 2.25,
            "toxicity": 0.7,
            "molecular_weight": 92.14,
        },
        {
            "smiles": "CC(C)O",
            "name": "isopropanol",
            "logP": 0.05,
            "solubility": -0.16,
            "toxicity": 0.0,
            "molecular_weight": 60.10,
        },
        {
            "smiles": "CCN(CC)CC",
            "name": "triethylamine",
            "logP": 1.45,
            "solubility": 0.30,
            "toxicity": 0.3,
            "molecular_weight": 101.19,
        },
        {
            "smiles": "CC(=O)OC",
            "name": "methyl_acetate",
            "logP": 0.18,
            "solubility": 0.00,
            "toxicity": 0.1,
            "molecular_weight": 74.08,
        },
        {
            "smiles": "CCCCCCCCCC",
            "name": "decane",
            "logP": 5.01,
            "solubility": 4.00,
            "toxicity": 0.2,
            "molecular_weight": 142.28,
        },
    ]
    
    # Validate SMILES
    validated_data = []
    for entry in sample_data:
        smiles = entry["smiles"]
        if validate_smiles(smiles):
            validated_data.append(entry)
        else:
            logger.warning(f"Invalid SMILES skipped: {smiles}")
    
    # Write to CSV
    if validated_data:
        fieldnames = ["smiles", "name", "logP", "solubility", "toxicity", "molecular_weight"]
        with open(output_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(validated_data)
        
        logger.info(f"Collected {len(validated_data)} molecules to {output_file}")
    else:
        logger.error("No valid molecules collected!")
    
    # TODO: Add ChEMBL API integration
    # TODO: Add ZINC database download
    # TODO: Add custom dataset loading


if __name__ == "__main__":
    collect_property_datasets()

