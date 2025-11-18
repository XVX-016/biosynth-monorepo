"""
SMILES featurization utilities using RDKit
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Optional, List


def featurize_smiles(smiles: str) -> Optional[List[float]]:
    """
    Convert SMILES string to Morgan fingerprint
    
    Args:
        smiles: SMILES string
        
    Returns:
        Feature vector (2048-bit fingerprint) or None if invalid SMILES
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    arr = list(fp.ToBitString())
    return [float(x) for x in arr]


def validate_smiles(smiles: str) -> bool:
    """
    Validate SMILES string
    
    Args:
        smiles: SMILES string to validate
        
    Returns:
        True if valid, False otherwise
    """
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None

