"""
RDKit Index - Fingerprint generation using existing featurizers

Reuses featurizers from Phase 5 (ml/featurize.py and ai/featurizer.py).
"""

from typing import Set, Optional
import logging

logger = logging.getLogger(__name__)


def compute_fingerprint(smiles: str, radius: int = 2, n_bits: int = 2048) -> Optional[Set[int]]:
    """
    Compute ECFP fingerprint for SMILES using existing featurizers.
    
    Tries to use RDKit via existing featurizers, falls back to simple hash-based.
    
    Args:
        smiles: SMILES string
        radius: ECFP radius (default: 2 for ECFP4)
        n_bits: Number of bits in fingerprint
    
    Returns:
        Set of active bit indices, or None if invalid SMILES
    """
    # Try RDKit-based featurizer first (from ai/featurizer.py)
    try:
        from backend.ai.featurizer import featurize_smiles
        
        fp_vector = featurize_smiles(smiles)
        if fp_vector is None:
            return None
        
        # Convert bit vector to set of active bit indices
        active_bits = {i for i, bit in enumerate(fp_vector) if bit > 0.5}
        return active_bits
    
    except ImportError:
        logger.warning("RDKit featurizer not available, trying alternative")
    
    # Try ml/featurize.py approach
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        # Generate ECFP fingerprint
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
        
        # Convert to set of active bit indices
        active_bits = set(fp.GetOnBits())
        return active_bits
    
    except ImportError:
        logger.warning("RDKit not available, using fallback fingerprint")
        return _fallback_fingerprint(smiles, n_bits)
    except Exception as e:
        logger.error(f"Error computing fingerprint: {e}")
        return None


def _fallback_fingerprint(smiles: str, n_bits: int) -> Set[int]:
    """
    Fallback fingerprint computation (simple hash-based).
    
    Used when RDKit is not available.
    """
    active_bits = set()
    
    # Hash different length substrings
    for length in [1, 2, 3, 4]:
        for i in range(len(smiles) - length + 1):
            substring = smiles[i:i+length]
            hash_val = hash(substring) % n_bits
            active_bits.add(hash_val)
    
    return active_bits


def validate_smiles(smiles: str) -> bool:
    """
    Validate SMILES string using existing validators.
    
    Args:
        smiles: SMILES string to validate
    
    Returns:
        True if valid, False otherwise
    """
    # Try existing validator
    try:
        from backend.ai.featurizer import validate_smiles
        return validate_smiles(smiles)
    except ImportError:
        pass
    
    # Fallback: try RDKit directly
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except ImportError:
        # If RDKit not available, basic check
        return len(smiles) > 0 and all(c.isalnum() or c in '()[]=+-.#@$%' for c in smiles)

