"""
SMARTS pattern parser and matcher
"""
from typing import Dict, List, Any, Optional
import re

class SMARTSParser:
    """Simple SMARTS parser for basic patterns"""
    
    def __init__(self):
        self.atom_classes = {
            'C': ['C'],
            'O': ['O'],
            'N': ['N'],
            'H': ['H'],
            'c': ['C'],  # aromatic carbon
            'o': ['O'],  # aromatic oxygen
            'n': ['N'],  # aromatic nitrogen
        }
    
    def parse(self, smarts: str) -> Dict[str, Any]:
        """
        Parse SMARTS string into pattern graph
        
        Returns:
            {
                "atoms": [...],
                "bonds": [...],
                "constraints": {...}
            }
        """
        # Simple parser for basic patterns like "C=O", "C-C", "c1ccccc1"
        atoms = []
        bonds = []
        constraints = {}
        
        # Tokenize SMARTS
        tokens = self._tokenize(smarts)
        
        # Build pattern (simplified - handles basic cases)
        atom_id = 0
        atom_map = {}
        
        i = 0
        while i < len(tokens):
            token = tokens[i]
            
            if token in self.atom_classes:
                # Create atom
                atom_id_str = f"a{atom_id}"
                atoms.append({
                    "id": atom_id_str,
                    "element": self.atom_classes[token][0],
                    "aromatic": token.islower()
                })
                atom_map[token] = atom_id_str
                atom_id += 1
            
            elif token in ['-', '=', '#', ':']:
                # Bond
                if i > 0 and i < len(tokens) - 1:
                    prev_atom = tokens[i-1]
                    next_atom = tokens[i+1]
                    if prev_atom in atom_map and next_atom in atom_map:
                        order = {'-': 1, '=': 2, '#': 3, ':': 1}[token]
                        bonds.append({
                            "atom1": atom_map[prev_atom],
                            "atom2": atom_map[next_atom],
                            "order": order
                        })
            
            i += 1
        
        return {
            "atoms": atoms,
            "bonds": bonds,
            "constraints": constraints
        }
    
    def _tokenize(self, smarts: str) -> List[str]:
        """Tokenize SMARTS string"""
        # Simple tokenizer - split by common patterns
        tokens = []
        i = 0
        while i < len(smarts):
            char = smarts[i]
            if char.isalnum():
                tokens.append(char)
            elif char in ['-', '=', '#', ':', '(', ')', '[', ']']:
                tokens.append(char)
            i += 1
        return tokens

def match_smarts(smarts_pattern: str, molecule: Dict[str, Any]) -> Tuple[bool, Optional[List[Dict[str, str]]]]:
    """
    Match SMARTS pattern against molecule
    
    Returns:
        (is_match, matched_atoms)
    """
    parser = SMARTSParser()
    pattern = parser.parse(smarts_pattern)
    
    # Convert pattern to molecule format for VF2
    pattern_mol = {
        "atoms": pattern["atoms"],
        "bonds": pattern["bonds"]
    }
    
    from .substructure import find_substructure
    is_match, mapping = find_substructure(pattern_mol, molecule)
    
    if is_match and mapping:
        matched_atoms = [{"pattern": p, "target": t} for p, t in mapping]
        return (True, matched_atoms)
    
    return (False, None)

