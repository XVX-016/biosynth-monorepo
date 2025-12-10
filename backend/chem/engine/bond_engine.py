"""
Bond Formation Engine - orchestrates ML prediction, rule validation, and optimization.
"""
from typing import Dict, List, Tuple
from ..core.models import Atom, Bond
from .rule_engine import RuleEngine
from .sanity_checker import MoleculeSanityChecker
from .optimizer import BondOptimizer
from ..ml.bond_ml_predictor import BondMLPredictor


class BondEngine:
    """
    Full bond formation pipeline:
    1. ML prediction for candidate bonds
    2. Rule-engine validation
    3. Sanity constraints (valency, charges)
    4. Optimization for most stable structure
    """
    
    def __init__(self):
        self.rule_engine = RuleEngine()
        self.ml_predictor = BondMLPredictor()
        self.sanity_checker = MoleculeSanityChecker()
        self.optimizer = BondOptimizer()

    def predict_all_bonds(
        self, 
        atoms: Dict[int, Atom], 
        positions: Dict[int, Tuple[float, float, float]]
    ) -> List[Bond]:
        """
        Full pipeline for bond prediction.
        
        Args:
            atoms: Dictionary of atoms by ID
            positions: Dictionary of 3D positions by atom ID
            
        Returns:
            List of optimized, validated bonds
        """
        # Step 1: ML prediction
        ml_bonds = self.ml_predictor.predict(atoms, positions)
        
        # Step 2: Rule validation
        valid_bonds = self.rule_engine.filter_invalid(atoms, ml_bonds)
        
        # Step 3: Sanity check (valence constraints)
        sane_bonds = self.sanity_checker.enforce(atoms, valid_bonds)
        
        # Step 4: Optimization
        optimized = self.optimizer.optimize(atoms, sane_bonds)

        return optimized

    def add_bond(self, atom_a: Atom, atom_b: Atom, order: int = 1) -> Bond:
        """
        Manual bond creation with validation.
        
        Args:
            atom_a: First atom
            atom_b: Second atom
            order: Bond order (default 1)
            
        Returns:
            Created bond
            
        Raises:
            ValueError: If bond violates chemistry rules
        """
        if not self.rule_engine.validate_bond(atom_a, atom_b, order):
            raise ValueError(f"Bond between {atom_a.element} and {atom_b.element} with order {order} violates standard rules.")
        
        return Bond(atom_a=atom_a.id, atom_b=atom_b.id, order=order)
