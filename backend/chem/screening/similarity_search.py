"""
Similarity Search Engine - High-level similarity search interface

Provides easy-to-use similarity search with multiple metrics and filters.
"""

from typing import List, Dict, Optional, Tuple
from .fingerprint_index import FingerprintIndex, get_fingerprint_index, compute_ecfp_fingerprint
import logging

logger = logging.getLogger(__name__)


class SimilaritySearchEngine:
    """
    High-level similarity search engine.
    
    Wraps fingerprint index with convenient search methods.
    """
    
    def __init__(self, index: Optional[FingerprintIndex] = None):
        """
        Args:
            index: FingerprintIndex instance (default: global instance)
        """
        self.index = index or get_fingerprint_index()
    
    def search_by_smiles(
        self,
        query_smiles: str,
        threshold: float = 0.5,
        max_results: int = 100,
        radius: int = 2
    ) -> List[Dict]:
        """
        Search for similar molecules by SMILES.
        
        Args:
            query_smiles: Query SMILES string
            threshold: Minimum similarity (0.0-1.0)
            max_results: Maximum results
            radius: ECFP radius
        
        Returns:
            List of result dicts with molecule_id, similarity, metadata
        """
        # Compute query fingerprint
        query_fp = compute_ecfp_fingerprint(query_smiles, radius=radius)
        
        if not query_fp:
            logger.warning(f"Failed to compute fingerprint for SMILES: {query_smiles}")
            return []
        
        # Search index
        results = self.index.search_similar(
            query_fp,
            threshold=threshold,
            max_results=max_results
        )
        
        # Format results with metadata
        formatted_results = []
        for molecule_id, similarity in results:
            metadata = self.index.get_metadata(molecule_id) or {}
            formatted_results.append({
                'molecule_id': molecule_id,
                'similarity': similarity,
                'smiles': metadata.get('smiles', ''),
                'name': metadata.get('name', ''),
                'metadata': metadata,
            })
        
        return formatted_results
    
    def search_by_fingerprint(
        self,
        query_fingerprint: set,
        threshold: float = 0.5,
        max_results: int = 100
    ) -> List[Dict]:
        """
        Search using pre-computed fingerprint.
        
        Args:
            query_fingerprint: Set of hash values
            threshold: Minimum similarity
            max_results: Maximum results
        
        Returns:
            List of result dicts
        """
        results = self.index.search_similar(
            query_fingerprint,
            threshold=threshold,
            max_results=max_results
        )
        
        formatted_results = []
        for molecule_id, similarity in results:
            metadata = self.index.get_metadata(molecule_id) or {}
            formatted_results.append({
                'molecule_id': molecule_id,
                'similarity': similarity,
                'smiles': metadata.get('smiles', ''),
                'name': metadata.get('name', ''),
                'metadata': metadata,
            })
        
        return formatted_results
    
    def get_k_nearest(
        self,
        query_smiles: str,
        k: int = 10,
        radius: int = 2
    ) -> List[Dict]:
        """
        Get k nearest neighbors.
        
        Args:
            query_smiles: Query SMILES
            k: Number of neighbors
            radius: ECFP radius
        
        Returns:
            List of k nearest neighbors
        """
        # Use very low threshold to get all candidates, then sort
        results = self.search_by_smiles(
            query_smiles,
            threshold=0.0,
            max_results=k * 10,  # Get more candidates
            radius=radius
        )
        
        # Sort by similarity and take top k
        results.sort(key=lambda x: x['similarity'], reverse=True)
        return results[:k]
    
    def batch_search(
        self,
        query_smiles_list: List[str],
        threshold: float = 0.5,
        max_results_per_query: int = 10
    ) -> Dict[str, List[Dict]]:
        """
        Batch similarity search for multiple queries.
        
        Args:
            query_smiles_list: List of query SMILES
            threshold: Minimum similarity
            max_results_per_query: Max results per query
        
        Returns:
            Dict mapping query SMILES to result lists
        """
        results = {}
        
        for query_smiles in query_smiles_list:
            results[query_smiles] = self.search_by_smiles(
                query_smiles,
                threshold=threshold,
                max_results=max_results_per_query
            )
        
        return results

