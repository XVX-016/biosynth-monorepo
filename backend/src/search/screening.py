"""
Screening Pipeline - Predicate-based molecular screening

Accepts any callable predicate(smiles) -> bool for flexible screening.
"""

from typing import List, Dict, Callable, Optional, Any
from .fingerprint_index import FingerprintIndex
from .rdkit_index import validate_smiles
import logging

logger = logging.getLogger(__name__)


class ScreeningPipeline:
    """
    Screening pipeline that accepts callable predicates.
    
    Supports property-threshold screening and custom predicates.
    """
    
    def __init__(self, index: Optional[FingerprintIndex] = None):
        """
        Args:
            index: FingerprintIndex instance to screen
        """
        self.index = index or FingerprintIndex()
    
    def batch_screen(
        self,
        library: List[str],
        predicate_fn: Callable[[str], bool],
        max_results: int = 1000
    ) -> List[Dict]:
        """
        Screen library using predicate function.
        
        Args:
            library: List of SMILES strings or molecule IDs
            predicate_fn: Callable that takes SMILES and returns bool
            max_results: Maximum results to return
        
        Returns:
            List of result dicts with molecule_id, smiles, passed, metadata
        """
        results = []
        count = 0
        
        for item in library:
            if count >= max_results:
                break
            
            # Handle both SMILES strings and molecule IDs
            if item in self.index.fingerprints:
                # It's a molecule ID
                metadata = self.index.get_metadata(item)
                smiles = metadata.get('smiles', '')
                molecule_id = item
            else:
                # Assume it's a SMILES string
                smiles = item
                molecule_id = None
                metadata = {}
            
            if not smiles:
                continue
            
            # Apply predicate
            try:
                passed = predicate_fn(smiles)
                if passed:
                    if molecule_id is None:
                        # Generate ID from SMILES hash
                        molecule_id = f"mol_{hash(smiles) % 1000000}"
                    
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
    
    def screen_index(
        self,
        predicate_fn: Callable[[str], bool],
        max_results: int = 1000
    ) -> List[Dict]:
        """
        Screen all molecules in index using predicate.
        
        Args:
            predicate_fn: Callable that takes SMILES and returns bool
            max_results: Maximum results to return
        
        Returns:
            List of result dicts
        """
        # Get all molecule IDs from index
        library = list(self.index.fingerprints.keys())
        return self.batch_screen(library, predicate_fn, max_results)
    
    def create_property_threshold_predicate(
        self,
        property_name: str,
        min_value: Optional[float] = None,
        max_value: Optional[float] = None
    ) -> Callable[[str], bool]:
        """
        Create predicate for property threshold screening.
        
        Args:
            property_name: Name of property to check
            min_value: Minimum value (inclusive)
            max_value: Maximum value (inclusive)
        
        Returns:
            Predicate function
        """
        # TODO: Integrate with property prediction from Phase 5
        # For now, use simple heuristics based on SMILES
        
        def predicate(smiles: str) -> bool:
            # Placeholder: would use actual property prediction
            # For MVP, use simple SMILES-based heuristics
            
            if property_name == 'molecular_weight':
                # Simple MW estimate from atom counts
                from backend.chem.ml.features import extract_features
                try:
                    # Would need molecule dict, so for now use fallback
                    # Count atoms in SMILES
                    atom_count = sum(1 for c in smiles if c.isupper() and c.isalpha())
                    # Rough MW estimate (C=12, O=16, N=14, etc.)
                    mw = atom_count * 15  # Rough average
                    
                    if min_value is not None and mw < min_value:
                        return False
                    if max_value is not None and mw > max_value:
                        return False
                    return True
                except:
                    return False
            
            elif property_name == 'logp':
                # Simple LogP estimate
                # Count polar vs non-polar atoms
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
                logger.warning(f"Unknown property: {property_name}")
                return True
        
        return predicate
    
    def create_custom_predicate(
        self,
        predicate_config: Dict[str, Any]
    ) -> Callable[[str], bool]:
        """
        Create custom predicate from configuration.
        
        Args:
            predicate_config: Configuration dict
        
        Returns:
            Predicate function
        """
        predicate_type = predicate_config.get('type')
        
        if predicate_type == 'property_threshold':
            return self.create_property_threshold_predicate(
                property_name=predicate_config.get('property_name'),
                min_value=predicate_config.get('min_value'),
                max_value=predicate_config.get('max_value'),
            )
        elif predicate_type == 'smiles_pattern':
            pattern = predicate_config.get('pattern', '')
            mode = predicate_config.get('mode', 'contains')
            
            def predicate(smiles: str) -> bool:
                if mode == 'contains':
                    return pattern in smiles
                elif mode == 'starts_with':
                    return smiles.startswith(pattern)
                elif mode == 'ends_with':
                    return smiles.endswith(pattern)
                elif mode == 'exact':
                    return smiles == pattern
                return False
            
            return predicate
        else:
            # Default: accept all
            logger.warning(f"Unknown predicate type: {predicate_type}")
            return lambda smiles: True

