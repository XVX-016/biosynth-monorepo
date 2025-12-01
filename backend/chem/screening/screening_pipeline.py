"""
Screening Pipeline - Multi-stage molecular screening workflow

Combines similarity search, substructure filtering, and property filters.
"""

from typing import List, Dict, Optional, Callable, Any
from .similarity_search import SimilaritySearchEngine
from .fingerprint_index import get_fingerprint_index
import logging

logger = logging.getLogger(__name__)


class ScreeningPipeline:
    """
    Multi-stage screening pipeline for virtual screening.
    
    Supports chaining multiple filters and search methods.
    """
    
    def __init__(self, similarity_engine: Optional[SimilaritySearchEngine] = None):
        """
        Args:
            similarity_engine: SimilaritySearchEngine instance
        """
        self.similarity_engine = similarity_engine or SimilaritySearchEngine()
        self.filters: List[Callable[[Dict], bool]] = []
    
    def add_filter(self, filter_func: Callable[[Dict], bool]):
        """
        Add a filter function to the pipeline.
        
        Args:
            filter_func: Function that takes molecule dict and returns bool
        """
        self.filters.append(filter_func)
    
    def add_property_filter(
        self,
        property_name: str,
        min_value: Optional[float] = None,
        max_value: Optional[float] = None
    ):
        """
        Add property-based filter.
        
        Args:
            property_name: Name of property to filter
            min_value: Minimum value (inclusive)
            max_value: Maximum value (inclusive)
        """
        def filter_func(mol_dict: Dict) -> bool:
            value = mol_dict.get('metadata', {}).get(property_name)
            if value is None:
                return False
            
            try:
                value = float(value)
            except (ValueError, TypeError):
                return False
            
            if min_value is not None and value < min_value:
                return False
            
            if max_value is not None and value > max_value:
                return False
            
            return True
        
        self.add_filter(filter_func)
    
    def add_smiles_filter(self, pattern: str, mode: str = 'contains'):
        """
        Add SMILES pattern filter.
        
        Args:
            pattern: SMILES pattern to match
            mode: 'contains', 'exact', 'starts_with', 'ends_with'
        """
        def filter_func(mol_dict: Dict) -> bool:
            smiles = mol_dict.get('smiles', '')
            
            if mode == 'contains':
                return pattern in smiles
            elif mode == 'exact':
                return smiles == pattern
            elif mode == 'starts_with':
                return smiles.startswith(pattern)
            elif mode == 'ends_with':
                return smiles.endswith(pattern)
            else:
                return False
        
        self.add_filter(filter_func)
    
    def screen_by_similarity(
        self,
        query_smiles: str,
        similarity_threshold: float = 0.5,
        max_results: int = 1000,
        apply_filters: bool = True
    ) -> List[Dict]:
        """
        Screen library by similarity to query molecule.
        
        Args:
            query_smiles: Query SMILES
            similarity_threshold: Minimum similarity
            max_results: Maximum results before filtering
            apply_filters: Apply registered filters
        
        Returns:
            List of screened molecules
        """
        # Initial similarity search
        candidates = self.similarity_engine.search_by_smiles(
            query_smiles,
            threshold=similarity_threshold,
            max_results=max_results
        )
        
        if not apply_filters or not self.filters:
            return candidates
        
        # Apply filters
        filtered = []
        for candidate in candidates:
            if all(filter_func(candidate) for filter_func in self.filters):
                filtered.append(candidate)
        
        return filtered
    
    def screen_by_substructure(
        self,
        query_smiles: str,
        similarity_threshold: float = 0.3,
        max_results: int = 1000,
        apply_filters: bool = True
    ) -> List[Dict]:
        """
        Screen library for molecules containing query as substructure.
        
        Uses lower similarity threshold since we're looking for substructures.
        
        Args:
            query_smiles: Query SMILES (substructure)
            similarity_threshold: Minimum similarity
            max_results: Maximum results
            apply_filters: Apply registered filters
        
        Returns:
            List of screened molecules
        """
        # Use similarity search with lower threshold for substructure matching
        candidates = self.similarity_engine.search_by_smiles(
            query_smiles,
            threshold=similarity_threshold,
            max_results=max_results
        )
        
        if not apply_filters or not self.filters:
            return candidates
        
        # Apply filters
        filtered = []
        for candidate in candidates:
            if all(filter_func(candidate) for filter_func in self.filters):
                filtered.append(candidate)
        
        return filtered
    
    def screen_batch(
        self,
        query_smiles_list: List[str],
        similarity_threshold: float = 0.5,
        max_results_per_query: int = 100,
        apply_filters: bool = True
    ) -> Dict[str, List[Dict]]:
        """
        Batch screening for multiple queries.
        
        Args:
            query_smiles_list: List of query SMILES
            similarity_threshold: Minimum similarity
            max_results_per_query: Max results per query
            apply_filters: Apply registered filters
        
        Returns:
            Dict mapping query SMILES to result lists
        """
        results = {}
        
        for query_smiles in query_smiles_list:
            results[query_smiles] = self.screen_by_similarity(
                query_smiles,
                similarity_threshold=similarity_threshold,
                max_results=max_results_per_query,
                apply_filters=apply_filters
            )
        
        return results
    
    def clear_filters(self):
        """Clear all registered filters."""
        self.filters = []
    
    def get_pipeline_stats(self) -> Dict:
        """Get pipeline statistics."""
        return {
            'num_filters': len(self.filters),
            'index_stats': get_fingerprint_index().get_stats(),
        }

