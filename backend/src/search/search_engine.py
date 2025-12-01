"""
Search Engine - High-level search interface

Provides similarity_search() and substructure_search() methods.
"""

from typing import List, Dict, Optional
from .fingerprint_index import FingerprintIndex
from .rdkit_index import compute_ecfp
from backend.chem.utils.validators import validate_smiles
from backend.chem.search.smarts import match_smarts
import logging

logger = logging.getLogger(__name__)


class SearchEngine:
    """
    High-level search engine for molecular search.
    """
    
    def __init__(self, index: Optional[FingerprintIndex] = None):
        """
        Args:
            index: FingerprintIndex instance (creates new if None)
        """
        from .fingerprint_index import FingerprintIndex
        self.index = index or FingerprintIndex()
    
    def similarity_search(
        self,
        smiles: str,
        k: int = 10,
        threshold: Optional[float] = None
    ) -> List[Dict]:
        """
        Search for similar molecules by SMILES.
        
        Args:
            smiles: Query SMILES string
            k: Number of results to return
            threshold: Minimum similarity threshold (optional)
        
        Returns:
            List of result dicts with molecule_id, smiles, similarity, metadata
        """
        # Validate SMILES
        if not validate_smiles(smiles):
            logger.warning(f"Invalid SMILES: {smiles}")
            return []
        
        # Compute query fingerprint
        query_fp = compute_ecfp(smiles)
        if not query_fp:
            logger.warning(f"Failed to compute fingerprint for: {smiles}")
            return []
        
        # Search index
        results = self.index.query_tanimoto(query_fp, k=k, threshold=threshold)
        
        # Format results
        formatted_results = []
        for molecule_id, similarity in results:
            metadata = self.index.get_metadata(molecule_id) or {}
            formatted_results.append({
                'molecule_id': molecule_id,
                'smiles': metadata.get('smiles', molecule_id),
                'similarity': similarity,
                'metadata': metadata,
            })
        
        return formatted_results
    
    def substructure_search(
        self,
        smarts: str,
        max_results: int = 100
    ) -> List[Dict]:
        """
        Search for molecules containing SMARTS pattern.
        
        Args:
            smarts: SMARTS pattern
            max_results: Maximum number of results
        
        Returns:
            List of result dicts with molecule_id, smiles, matched_atoms, metadata
        """
        results = []
        count = 0
        
        # Search all molecules in index
        for molecule_id in self.index.fingerprints.keys():
            if count >= max_results:
                break
            
            metadata = self.index.get_metadata(molecule_id)
            smiles = metadata.get('smiles', '')
            
            if not smiles:
                continue
            
            # Convert SMILES to molecule dict format for match_smarts
            # For now, we'll use a simple approach
            # TODO: Properly convert SMILES to molecule dict format
            try:
                # Use existing SMARTS matcher
                # Note: match_smarts expects molecule dict with atoms/bonds
                # For MVP, we'll do a simple text-based check
                # In production, would convert SMILES to molecule dict
                
                # Simple fallback: check if SMARTS pattern appears in SMILES
                # This is not correct SMARTS matching, but works for MVP
                if smarts in smiles:
                    results.append({
                        'molecule_id': molecule_id,
                        'smiles': smiles,
                        'matched_atoms': [],  # Would need proper matching
                        'metadata': metadata,
                    })
                    count += 1
            except Exception as e:
                logger.warning(f"Error matching SMARTS for {molecule_id}: {e}")
                continue
        
        return results
    
    def get_index(self) -> FingerprintIndex:
        """Get the underlying fingerprint index."""
        return self.index
