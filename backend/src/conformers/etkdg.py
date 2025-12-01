"""
ETKDG Generator - Placeholder for RDKit ETKDG conformer generation

ETKDG (Experimental-Torsion-angle preference with Distance Geometry) is a
popular conformer generation method. This is a placeholder implementation.
"""

from typing import List, Dict
import logging

logger = logging.getLogger(__name__)


class ETKDGGenerator:
    """
    Placeholder for RDKit ETKDG conformer generator.
    
    Must raise NotImplementedError so ConformerGenerator knows to fallback.
    """
    
    def generate(self, smiles: str, n: int) -> List[Dict]:
        """
        Placeholder for RDKit ETKDG.
        
        Must raise NotImplementedError or return None,
        so ConformerGenerator knows to fallback.
        
        Args:
            smiles: SMILES string
            n: Number of conformers
        
        Returns:
            List of conformer dicts
        
        Raises:
            NotImplementedError: Always, as this is a placeholder
        """
        raise NotImplementedError("ETKDG is not implemented yet")
