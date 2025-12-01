"""
Fingerprint Index - Fast molecular fingerprint storage and retrieval

Uses ECFP (Extended Connectivity Fingerprints) for similarity search.
"""

import hashlib
from typing import Dict, List, Set, Optional, Tuple
from collections import defaultdict
import logging

logger = logging.getLogger(__name__)


class FingerprintIndex:
    """
    Index for molecular fingerprints (ECFP-based).
    
    Enables fast similarity search by storing fingerprints and
    maintaining inverted indices for common substructures.
    """
    
    def __init__(self):
        # molecule_id -> fingerprint (set of hash values)
        self.fingerprints: Dict[str, Set[int]] = {}
        
        # hash_value -> set of molecule_ids (inverted index)
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
        Add molecule fingerprint to index.
        
        Args:
            molecule_id: Unique molecule identifier
            fingerprint: Set of hash values (ECFP bits)
            metadata: Optional metadata (SMILES, name, etc.)
        """
        # Remove old fingerprint if exists
        if molecule_id in self.fingerprints:
            self.remove_molecule(molecule_id)
        
        # Add new fingerprint
        self.fingerprints[molecule_id] = fingerprint.copy()
        
        # Update inverted index
        for hash_val in fingerprint:
            self.inverted_index[hash_val].add(molecule_id)
        
        # Store metadata
        self.metadata[molecule_id] = metadata or {}
    
    def remove_molecule(self, molecule_id: str):
        """Remove molecule from index."""
        if molecule_id not in self.fingerprints:
            return
        
        # Remove from inverted index
        fingerprint = self.fingerprints[molecule_id]
        for hash_val in fingerprint:
            self.inverted_index[hash_val].discard(molecule_id)
            # Clean up empty sets
            if not self.inverted_index[hash_val]:
                del self.inverted_index[hash_val]
        
        # Remove from fingerprints and metadata
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
        Find candidate molecules that share at least min_common_bits with query.
        
        Uses inverted index for fast candidate filtering.
        """
        if not query_fingerprint:
            return set()
        
        # Count occurrences per molecule
        candidate_counts: Dict[str, int] = defaultdict(int)
        
        for hash_val in query_fingerprint:
            if hash_val in self.inverted_index:
                for molecule_id in self.inverted_index[hash_val]:
                    candidate_counts[molecule_id] += 1
        
        # Filter by minimum common bits
        candidates = {
            mol_id for mol_id, count in candidate_counts.items()
            if count >= min_common_bits
        }
        
        return candidates
    
    def search_similar(
        self,
        query_fingerprint: Set[int],
        threshold: float = 0.5,
        max_results: int = 100
    ) -> List[Tuple[str, float]]:
        """
        Search for similar molecules using Tanimoto similarity.
        
        Args:
            query_fingerprint: Query molecule fingerprint
            threshold: Minimum similarity score (0.0-1.0)
            max_results: Maximum number of results
        
        Returns:
            List of (molecule_id, similarity_score) tuples, sorted by score
        """
        # Find candidates using inverted index
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
            
            if similarity >= threshold:
                results.append((molecule_id, similarity))
        
        # Sort by similarity (descending)
        results.sort(key=lambda x: x[1], reverse=True)
        
        # Limit results
        return results[:max_results]
    
    def get_stats(self) -> Dict:
        """Get index statistics."""
        return {
            'num_molecules': len(self.fingerprints),
            'num_fingerprint_bits': len(self.inverted_index),
            'avg_fingerprint_size': (
                sum(len(fp) for fp in self.fingerprints.values()) / len(self.fingerprints)
                if self.fingerprints else 0
            ),
        }


# Global fingerprint index instance
_fingerprint_index: Optional[FingerprintIndex] = None


def get_fingerprint_index() -> FingerprintIndex:
    """Get global fingerprint index instance."""
    global _fingerprint_index
    if _fingerprint_index is None:
        _fingerprint_index = FingerprintIndex()
    return _fingerprint_index


def compute_ecfp_fingerprint(
    smiles: str,
    radius: int = 2,
    n_bits: int = 2048
) -> Set[int]:
    """
    Compute ECFP fingerprint for SMILES string.
    
    This is a simplified implementation. In production, use RDKit.
    
    Args:
        smiles: SMILES string
        radius: ECFP radius (default: 2 for ECFP4)
        n_bits: Number of bits in fingerprint
    
    Returns:
        Set of hash values (active bits)
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return set()
        
        # Generate ECFP fingerprint
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
        
        # Convert to set of active bit indices
        active_bits = set(fp.GetOnBits())
        return active_bits
    
    except ImportError:
        # Fallback: simple hash-based fingerprint
        logger.warning("RDKit not available, using fallback fingerprint")
        return _fallback_fingerprint(smiles, n_bits)


def _fallback_fingerprint(smiles: str, n_bits: int) -> Set[int]:
    """
    Fallback fingerprint computation (simple hash-based).
    
    This is a mock implementation for when RDKit is not available.
    """
    # Simple approach: hash substrings
    active_bits = set()
    
    # Hash different length substrings
    for length in [1, 2, 3]:
        for i in range(len(smiles) - length + 1):
            substring = smiles[i:i+length]
            hash_val = hash(substring) % n_bits
            active_bits.add(hash_val)
    
    return active_bits

