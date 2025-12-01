"""
SMILES validation utilities
"""

from backend.ai.featurizer import validate_smiles as _validate_smiles

# Re-export for convenience
def validate_smiles(smiles: str) -> bool:
    """
    Validate SMILES string.
    
    Args:
        smiles: SMILES string to validate
    
    Returns:
        True if valid, False otherwise
    """
    return _validate_smiles(smiles)

