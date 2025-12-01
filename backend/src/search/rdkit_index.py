"""
RDKit Index - ECFP fingerprint computation

Uses existing featurizer to compute fingerprints.
"""

from typing import Optional, Set
from backend.ai.featurizer import get_ecfp


def compute_ecfp(smiles: str, radius: int = 2, n_bits: int = 2048) -> Optional[Set[int]]:
    """
    Compute ECFP fingerprint for SMILES.
    
    Args:
        smiles: SMILES string
        radius: ECFP radius (default: 2 for ECFP4)
        n_bits: Number of bits in fingerprint
    
    Returns:
        Set of active bit indices, or None if invalid SMILES
    """
    return get_ecfp(smiles, radius=radius, n_bits=n_bits)
