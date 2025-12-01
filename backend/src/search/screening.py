"""
Screening Pipeline - Generic predicate-based screening

Accepts any callable predicate(smiles) -> bool for flexible screening.
"""

from typing import List, Dict, Callable, Optional
from .fingerprint_index import FingerprintIndex
import logging

logger = logging.getLogger(__name__)


class ScreeningPipeline:
    """
    Screening pipeline that accepts callable predicates.
    """
    
    def __init__(self, index: Optional[FingerprintIndex] = None):
        """
        Args:
            index: FingerprintIndex instance to screen
        """
        from .fingerprint_index import FingerprintIndex
        self.index = index or FingerprintIndex()
    
    def property_threshold_predicate(
        self,
        name: str,
        min_value: Optional[float] = None,
        max_value: Optional[float] = None
    ) -> Callable[[str], bool]:
        """
        Create predicate for property threshold screening.
        
        Args:
            name: Property name (e.g., 'logp', 'molecular_weight')
            min_value: Minimum value (inclusive)
            max_value: Maximum value (inclusive)
        
        Returns:
            Predicate function that takes SMILES and returns bool
        """
        def predicate(smiles: str) -> bool:
            # TODO: Integrate with property prediction from Phase 5
            # For now, use simple heuristics based on SMILES
            
            if name == 'molecular_weight':
                # Simple MW estimate from atom counts
                atom_count = sum(1 for c in smiles if c.isupper() and c.isalpha())
                mw = atom_count * 15  # Rough average
            
                if min_value is not None and mw < min_value:
                    return False
                if max_value is not None and mw > max_value:
                    return False
                return True
            
            elif name == 'logp':
                # Simple LogP estimate
                polar_atoms = sum(1 for c in smiles if c in 'ONFClBrI')
                non_polar_atoms = sum(1 for c in smiles if c == 'C')
                logp = 0.5 * non_polar_atoms - 0.5 * polar_atoms
            
                if min_value is not None and logp < min_value:
                    return False
                if max_value is not None and logp > max_value:
                    return False
                return True
            
            else:
                # Unknown property, return True (pass all)
                logger.warning(f"Unknown property: {name}")
                return True
        
        return predicate
    
    def batch_screen(
        self,
        pred: Callable[[str], bool],
        max_results: int = 1000
    ) -> List[Dict]:
        """
        Screen library using predicate function.
        
        Args:
            pred: Callable that takes SMILES and returns bool
            max_results: Maximum results to return
        
        Returns:
            List of result dicts with molecule_id, smiles, passed, metadata
        """
        results = []
        count = 0
        
        for molecule_id in self.index.fingerprints.keys():
            if count >= max_results:
                break
            
            metadata = self.index.get_metadata(molecule_id)
            smiles = metadata.get('smiles', '')
            
            if not smiles:
                continue
            
            # Apply predicate
            try:
                passed = pred(smiles)
                if passed:
                    results.append({
                        'molecule_id': molecule_id,
                        'smiles': smiles,
                        'passed': True,
                        'metadata': metadata,
                    })
                    count += 1
            except Exception as e:
                logger.warning(f"Error applying predicate to {smiles}: {e}")
                continue
        
        return results
    
    def get_index(self) -> FingerprintIndex:
        """Get the fingerprint index."""
        return self.index
