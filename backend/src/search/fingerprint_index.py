"""
Fingerprint Index - Simple in-memory fingerprint storage

Reuses existing featurizers from Phase 5 for fingerprint generation.
"""

from typing import Dict, Set, List, Optional, Tuple
from collections import defaultdict
import logging

logger = logging.getLogger(__name__)


class FingerprintIndex:
    """
    Simple in-memory fingerprint index for similarity search.
    
    Uses ECFP fingerprints computed via existing featurizers.
    """
    
    def __init__(self):
        # molecule_id -> fingerprint (set of active bits)
        self.fingerprints: Dict[str, Set[int]] = {}
        
        # bit_index -> set of molecule_ids (inverted index for fast lookup)
        self.inverted_index: Dict[int, Set[str]] = defaultdict(set)
        
        # molecule_id -> metadata
        self.metadata: Dict[str, Dict] = {}
    
    def add_molecule(
        self,
        molecule_id: str,
        fingerprint: Set[int],
        metadata: Optional[Dict] = None
    ):
        """
        Add molecule to index.
        
        Args:
            molecule_id: Unique identifier
            fingerprint: Set of active bit indices
            metadata: Optional metadata dict
        """
        # Remove old entry if exists
        if molecule_id in self.fingerprints:
            self.remove_molecule(molecule_id)
        
        # Add fingerprint
        self.fingerprints[molecule_id] = fingerprint.copy()
        
        # Update inverted index
        for bit_index in fingerprint:
            self.inverted_index[bit_index].add(molecule_id)
        
        # Store metadata
        self.metadata[molecule_id] = metadata or {}
    
    def remove_molecule(self, molecule_id: str):
        """Remove molecule from index."""
        if molecule_id not in self.fingerprints:
            return
        
        # Remove from inverted index
        fingerprint = self.fingerprints[molecule_id]
        for bit_index in fingerprint:
            self.inverted_index[bit_index].discard(molecule_id)
            # Clean up empty sets
            if not self.inverted_index[bit_index]:
                del self.inverted_index[bit_index]
        
        # Remove from storage
        del self.fingerprints[molecule_id]
        if molecule_id in self.metadata:
            del self.metadata[molecule_id]
    
    def get_fingerprint(self, molecule_id: str) -> Optional[Set[int]]:
        """Get fingerprint for molecule."""
        return self.fingerprints.get(molecule_id)
    
    def get_metadata(self, molecule_id: str) -> Optional[Dict]:
        """Get metadata for molecule."""
        return self.metadata.get(molecule_id)
    
    def compute_tanimoto_similarity(
        self,
        fp1: Set[int],
        fp2: Set[int]
    ) -> float:
        """
        Compute Tanimoto coefficient (Jaccard similarity).
        
        Returns:
            Similarity score between 0.0 and 1.0
        """
        if not fp1 or not fp2:
            return 0.0
        
        intersection = len(fp1 & fp2)
        union = len(fp1 | fp2)
        
        if union == 0:
            return 0.0
        
        return intersection / union
    
    def find_candidates(
        self,
        query_fingerprint: Set[int],
        min_common_bits: int = 1
    ) -> Set[str]:
        """
        Find candidate molecules sharing at least min_common_bits.
        
        Uses inverted index for fast filtering.
        """
        if not query_fingerprint:
            return set()
        
        # Count common bits per molecule
        candidate_counts: Dict[str, int] = defaultdict(int)
        
        for bit_index in query_fingerprint:
            if bit_index in self.inverted_index:
                for molecule_id in self.inverted_index[bit_index]:
                    candidate_counts[molecule_id] += 1
        
        # Filter by minimum common bits
        return {
            mol_id for mol_id, count in candidate_counts.items()
            if count >= min_common_bits
        }
    
    def search_similar(
        self,
        query_fingerprint: Set[int],
        k: int = 10,
        threshold: Optional[float] = None
    ) -> List[Tuple[str, float]]:
        """
        Search for similar molecules.
        
        Args:
            query_fingerprint: Query fingerprint (set of bit indices)
            k: Number of results
            threshold: Minimum similarity (optional)
        
        Returns:
            List of (molecule_id, similarity) tuples, sorted by similarity
        """
        # Find candidates
        candidates = self.find_candidates(query_fingerprint)
        
        if not candidates:
            return []
        
        # Compute similarities
        results = []
        for molecule_id in candidates:
            target_fp = self.fingerprints.get(molecule_id)
            if not target_fp:
                continue
            
            similarity = self.compute_tanimoto_similarity(
                query_fingerprint,
                target_fp
            )
            
            # Apply threshold if specified
            if threshold is not None and similarity < threshold:
                continue
            
            results.append((molecule_id, similarity))
        
        # Sort by similarity (descending)
        results.sort(key=lambda x: x[1], reverse=True)
        
        # Return top k
        return results[:k]
    
    def get_stats(self) -> Dict:
        """Get index statistics."""
        return {
            'num_molecules': len(self.fingerprints),
            'num_bits': len(self.inverted_index),
            'avg_fingerprint_size': (
                sum(len(fp) for fp in self.fingerprints.values()) / len(self.fingerprints)
                if self.fingerprints else 0
            ),
        }
    
    def clear(self):
        """Clear all molecules from index."""
        self.fingerprints.clear()
        self.inverted_index.clear()
        self.metadata.clear()

