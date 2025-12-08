"""
ML-based bond prediction wrapper.
"""
from typing import Dict, List, Tuple
from ..core.models import Atom, Bond
from .bond_predictor import BondPredictor


class BondMLPredictor:
    """Wrapper around ML bond predictor for engine integration."""
    
    def __init__(self):
        self.predictor = BondPredictor()
    
    def predict(self, atoms: Dict[int, Atom], positions: Dict[int, Tuple[float, float, float]]) -> List[Bond]:
        """
        Predict bonds using ML model or heuristic fallback.
        
        Args:
            atoms: Dictionary of atoms by ID
            positions: Dictionary of positions by atom ID
            
        Returns:
            List of predicted bonds
        """
        # Convert to format expected by BondPredictor
        atom_list = []
        for atom_id, atom in atoms.items():
            pos = positions.get(atom_id, (0.0, 0.0, 0.0))
            atom_list.append((str(atom_id), atom.element, pos))
        
        # Get predictions
        bond_dicts = self.predictor.predict(atom_list)
        
        # Convert to Bond objects
        bonds = []
        for b in bond_dicts:
            try:
                bonds.append(Bond(
                    atom_a=int(b['a']),
                    atom_b=int(b['b']),
                    order=b.get('order', 1)
                ))
            except (KeyError, ValueError) as e:
                # Skip malformed bonds
                continue
        
        return bonds
