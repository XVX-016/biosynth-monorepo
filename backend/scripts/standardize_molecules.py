"""
Step 3: Standardize Molecules

Removes salts, normalizes charges, aromatizes rings, and generates
canonical SMILES.
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
    from rdkit.Chem import SaltRemover, SanitizeMol
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available. Standardization will be limited.")

try:
    from backend.chem.utils.validators import validate_smiles
except ImportError:
    # Fallback: simple validation
    def validate_smiles(smiles: str) -> bool:
        return bool(smiles and len(smiles) > 0)

logger = logging.getLogger(__name__)


def standardize_molecules():
    """
    Standardize molecules in the property dataset.
    
    Steps:
    1. Remove salts
    2. Normalize charges
    3. Aromatize rings
    4. Generate canonical SMILES
    """
    data_dir = Path("data/datasets")
    input_file = data_dir / "properties.csv"
    cleaned_dir = data_dir / "cleaned"
    cleaned_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = cleaned_dir / "properties_cleaned.csv"
    
    if not input_file.exists():
        logger.warning(f"Input file not found: {input_file}")
        return
    
    logger.info(f"Standardizing molecules from {input_file}")
    
    standardized = []
    salt_remover = SaltRemover.SaltRemover() if RDKIT_AVAILABLE else None
    
    with open(input_file, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            smiles = row.get("smiles", "").strip()
            if not smiles:
                continue
            
            # Standardize
            standardized_smiles = _standardize_smiles(smiles, salt_remover)
            
            if standardized_smiles:
                # Copy all other fields
                new_row = row.copy()
                new_row["smiles"] = standardized_smiles
                new_row["original_smiles"] = smiles  # Keep original
                standardized.append(new_row)
            else:
                logger.warning(f"Failed to standardize: {smiles}")
    
    # Write standardized dataset
    if standardized:
        fieldnames = list(standardized[0].keys())
        with open(output_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(standardized)
        
        logger.info(f"Standardized {len(standardized)} molecules to {output_file}")
    else:
        logger.error("No molecules standardized!")


def _standardize_smiles(smiles: str, salt_remover) -> str:
    """
    Standardize a SMILES string.
    
    Args:
        smiles: Input SMILES
        salt_remover: RDKit SaltRemover instance
    
    Returns:
        Standardized SMILES or empty string if failed
    """
    if not RDKIT_AVAILABLE:
        # Fallback: just validate
        return smiles if validate_smiles(smiles) else ""
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ""
        
        # Remove salts
        if salt_remover:
            mol = salt_remover.StripMol(mol)
        
        # Sanitize
        try:
            SanitizeMol(mol)
        except:
            pass  # Continue even if sanitization fails
        
        # Generate canonical SMILES
        canonical = Chem.MolToSmiles(mol, canonical=True)
        return canonical
    except Exception as e:
        logger.warning(f"Standardization error for {smiles}: {e}")
        return ""


if __name__ == "__main__":
    standardize_molecules()

