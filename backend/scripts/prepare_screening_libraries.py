"""
Step 2: Prepare Screening Libraries

Extracts subsets from property dataset and saves as .smi files
for Phase 7 screening.
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
    from backend.chem.utils.validators import validate_smiles
except ImportError:
    # Fallback: simple validation
    def validate_smiles(smiles: str) -> bool:
        return bool(smiles and len(smiles) > 0)

logger = logging.getLogger(__name__)


def prepare_screening_libraries():
    """
    Prepare screening libraries from property datasets.
    
    Creates .smi files in /data/libraries/ with format:
    SMILES <space> molecule_name
    """
    data_dir = Path("data")
    libraries_dir = data_dir / "libraries"
    libraries_dir.mkdir(parents=True, exist_ok=True)
    
    properties_file = data_dir / "datasets" / "properties.csv"
    
    if not properties_file.exists():
        logger.warning(f"Properties file not found: {properties_file}")
        logger.info("Creating sample library from scratch...")
        _create_sample_library(libraries_dir)
        return
    
    logger.info(f"Reading properties from {properties_file}")
    
    # Read properties CSV
    molecules = []
    with open(properties_file, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            smiles = row.get("smiles", "").strip()
            name = row.get("name", "").strip() or smiles
            
            if validate_smiles(smiles):
                molecules.append((smiles, name))
            else:
                logger.warning(f"Invalid SMILES skipped: {smiles}")
    
    if not molecules:
        logger.warning("No valid molecules found. Creating sample library...")
        _create_sample_library(libraries_dir)
        return
    
    # Write to .smi file
    output_file = libraries_dir / "compounds.smi"
    with open(output_file, "w") as f:
        for smiles, name in molecules:
            f.write(f"{smiles} {name}\n")
    
    logger.info(f"Created screening library: {output_file} ({len(molecules)} molecules)")
    
    # Optionally create smaller subsets
    if len(molecules) > 100:
        # Create small subset for testing
        subset_file = libraries_dir / "compounds_small.smi"
        with open(subset_file, "w") as f:
            for smiles, name in molecules[:50]:
                f.write(f"{smiles} {name}\n")
        logger.info(f"Created small subset: {subset_file} (50 molecules)")


def _create_sample_library(libraries_dir: Path):
    """Create a sample library if no properties file exists."""
    sample_molecules = [
        ("CCO", "ethanol"),
        ("CC(=O)O", "acetic_acid"),
        ("c1ccccc1", "benzene"),
        ("CCc1ccccc1", "toluene"),
        ("CC(C)O", "isopropanol"),
        ("CCN(CC)CC", "triethylamine"),
        ("CC(=O)OC", "methyl_acetate"),
        ("CCCCCCCCCC", "decane"),
    ]
    
    output_file = libraries_dir / "compounds.smi"
    with open(output_file, "w") as f:
        for smiles, name in sample_molecules:
            f.write(f"{smiles} {name}\n")
    
    logger.info(f"Created sample library: {output_file} ({len(sample_molecules)} molecules)")


if __name__ == "__main__":
    prepare_screening_libraries()

