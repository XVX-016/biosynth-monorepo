"""
Fingerprint Index - Simple in-memory fingerprint storage

Uses simple Python lists for storage.
"""

from typing import List, Dict, Set, Optional, Tuple


class FingerprintIndex:
    """
    Simple in-memory fingerprint index.
    
    Stores fingerprints and metadata for fast similarity search.
    """
    
    def __init__(self):
        # molecule_id -> fingerprint (set of active bits)
        self.fingerprints: Dict[str, Set[int]] = {}
        # molecule_id -> metadata
        self.metadata: Dict[str, Dict] = {}
    
    def add(self, smiles: str, fp: Set[int], metadata: Optional[Dict] = None):
        """
        Add molecule to index.
        
        Args:
            smiles: SMILES string (used as ID if no ID in metadata)
            fp: Fingerprint (set of active bit indices)
            metadata: Optional metadata dict
        """
        # Use smiles as ID if no ID provided
        molecule_id = metadata.get('id', smiles) if metadata else smiles
        
        self.fingerprints[molecule_id] = fp.copy()
        self.metadata[molecule_id] = metadata or {'smiles': smiles}
        if 'smiles' not in self.metadata[molecule_id]:
            self.metadata[molecule_id]['smiles'] = smiles
    
    def query_tanimoto(
        self,
        fp: Set[int],
        k: int = 10,
        threshold: Optional[float] = None
    ) -> List[Tuple[str, float]]:
        """
        Query index using Tanimoto similarity.
        
        Args:
            fp: Query fingerprint (set of active bits)
            k: Number of results to return
            threshold: Minimum similarity threshold (optional)
        
        Returns:
            List of (molecule_id, similarity) tuples, sorted by similarity
        """
        results = []
        
        for molecule_id, target_fp in self.fingerprints.items():
            # Compute Tanimoto coefficient
            intersection = len(fp & target_fp)
            union = len(fp | target_fp)
            
            if union == 0:
                similarity = 0.0
            else:
                similarity = intersection / union
            
            # Apply threshold if specified
            if threshold is not None and similarity < threshold:
                continue
            
            results.append((molecule_id, similarity))
        
        # Sort by similarity (descending)
        results.sort(key=lambda x: x[1], reverse=True)
        
        # Return top k
        return results[:k]
    
    def get_metadata(self, molecule_id: str) -> Optional[Dict]:
        """Get metadata for molecule."""
        return self.metadata.get(molecule_id)
    
    def clear(self):
        """Clear all molecules from index."""
        self.fingerprints.clear()
        self.metadata.clear()
    
    def size(self) -> int:
        """Get number of molecules in index."""
        return len(self.fingerprints)
