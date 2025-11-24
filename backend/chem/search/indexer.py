"""
Molecule indexer for fast search
"""
from typing import Dict, List, Any
from .tokenizer import tokenize_molecule

class MoleculeIndex:
    """Index for fast substructure search"""
    
    def __init__(self):
        self.index: Dict[str, List[str]] = {}  # fingerprint -> molecule_ids
        self.molecules: Dict[str, Dict[str, Any]] = {}  # molecule_id -> tokenized
    
    def add_molecule(self, molecule_id: str, molecule: Dict[str, Any]):
        """Add molecule to index"""
        tokenized = tokenize_molecule(molecule)
        self.molecules[molecule_id] = tokenized
        
        # Index by fingerprint
        fp = tokenized["fingerprint"]
        if fp not in self.index:
            self.index[fp] = []
        self.index[fp].append(molecule_id)
    
    def search_candidates(self, query_molecule: Dict[str, Any]) -> List[str]:
        """
        Find candidate molecules that might match query
        Returns list of molecule IDs to check with VF2
        """
        query_tokenized = tokenize_molecule(query_molecule)
        query_fp = query_tokenized["fingerprint"]
        
        # Simple candidate selection: molecules with similar fingerprints
        candidates = []
        for fp, mol_ids in self.index.items():
            # If query fingerprint is substring or similar, add as candidate
            if query_fp in fp or fp in query_fp:
                candidates.extend(mol_ids)
        
        # If no candidates, return all (fallback)
        if not candidates:
            candidates = list(self.molecules.keys())
        
        return candidates
    
    def get_molecule(self, molecule_id: str) -> Dict[str, Any]:
        """Get tokenized molecule by ID"""
        return self.molecules.get(molecule_id, {})

# Global index instance
_global_index = MoleculeIndex()

def get_index() -> MoleculeIndex:
    """Get global molecule index"""
    return _global_index

