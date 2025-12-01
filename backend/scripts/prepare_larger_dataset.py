"""
Prepare Larger Dataset for Training

Generates a larger synthetic dataset for better model training.
"""

import sys
from pathlib import Path
import csv
import random
import logging

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

from rdkit import Chem
from rdkit.Chem import Descriptors
import numpy as np

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def generate_synthetic_dataset(n_molecules=1000):
    """Generate synthetic molecular dataset."""
    logger.info(f"Generating {n_molecules} synthetic molecules...")
    
    # Common molecular fragments
    fragments = [
        "C", "CC", "CCC", "CCCC", "CCO", "CCCO", "CC(C)O",
        "c1ccccc1", "CCc1ccccc1", "c1ccccc1O", "c1ccccc1C",
        "CC(=O)O", "CC(=O)C", "CCN", "CCCN", "CC(C)N",
        "CCOC", "CCCOC", "CC(C)OC", "CCOCC",
    ]
    
    molecules = []
    seen_smiles = set()
    
    # Generate combinations
    for _ in range(n_molecules):
        # Random combination of fragments
        n_fragments = random.randint(1, 3)
        selected = random.sample(fragments, min(n_fragments, len(fragments)))
        
        # Try to combine (simplified - just concatenate)
        smiles = "".join(selected)
        
        # Validate with RDKit
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol and smiles not in seen_smiles:
                canonical = Chem.MolToSmiles(mol)
                if canonical not in seen_smiles:
                    seen_smiles.add(canonical)
                    molecules.append(canonical)
        except:
            continue
        
        if len(molecules) >= n_molecules:
            break
    
    logger.info(f"Generated {len(molecules)} valid unique molecules")
    return molecules


def compute_properties(smiles_list):
    """Compute molecular properties."""
    logger.info("Computing properties...")
    
    data = []
    for smiles in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                continue
            
            # Compute properties
            logP = Descriptors.MolLogP(mol)
            mw = Descriptors.MolWt(mol)
            solubility = -np.log10(max(0.001, Descriptors.MolMR(mol) / 100.0))  # Mock solubility
            toxicity = random.uniform(0.0, 1.0)  # Mock toxicity (would use real model)
            
            data.append({
                "smiles": smiles,
                "logP": logP,
                "molecular_weight": mw,
                "solubility": solubility,
                "toxicity": toxicity,
            })
        except Exception as e:
            logger.warning(f"Failed to compute properties for {smiles}: {e}")
            continue
    
    logger.info(f"Computed properties for {len(data)} molecules")
    return data


def save_dataset(data, output_file):
    """Save dataset to CSV."""
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["smiles", "logP", "molecular_weight", "solubility", "toxicity"])
        writer.writeheader()
        writer.writerows(data)
    
    logger.info(f"Saved dataset to {output_path}")


def main():
    """Main function."""
    logger.info("=" * 60)
    logger.info("Prepare Larger Dataset")
    logger.info("=" * 60)
    
    # Generate molecules
    molecules = generate_synthetic_dataset(n_molecules=1000)
    
    # Compute properties
    data = compute_properties(molecules)
    
    # Save
    save_dataset(data, "data/datasets/properties_large.csv")
    
    logger.info("\n" + "=" * 60)
    logger.info("Dataset preparation complete!")
    logger.info("=" * 60)
    logger.info(f"\nNext steps:")
    logger.info(f"1. Run: python scripts/prepare_datasets.py")
    logger.info(f"2. Run: python scripts/train_ml_models.py")


if __name__ == "__main__":
    main()

