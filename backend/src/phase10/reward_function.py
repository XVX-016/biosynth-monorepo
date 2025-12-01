"""
Reward Function - Composite scoring from ML predictions, screening, QM/MD

Combines multiple property predictions and constraints into
a single reward score for RL optimization.
"""

from typing import List, Dict, Optional, Any, Callable
import logging

logger = logging.getLogger(__name__)


class RewardFunction:
    """
    Composite reward function for molecule optimization.
    
    Combines:
    - ML property predictions (logP, solubility, toxicity, etc.)
    - Screening results (similarity, substructure matches)
    - QM/MD properties (energy, stability)
    - Custom constraints (molecular weight, drug-likeness, etc.)
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize reward function.
        
        Args:
            config: Configuration dict with:
                - weights: Dict mapping property names to weights
                - constraints: Dict of constraint functions
                - normalization: Whether to normalize scores
        """
        self.config = config or {
            "weights": {
                "logP": 0.2,
                "solubility": 0.2,
                "toxicity": -0.3,  # Negative (lower is better)
                "similarity": 0.1,
                "energy": -0.1,  # Negative (lower is better)
                "drug_likeness": 0.1,
            },
            "constraints": {},
            "normalization": True,
        }
        logger.info("Reward Function initialized")
    
    def compute(
        self,
        smiles: str,
        ml_predictions: Optional[Dict[str, float]] = None,
        screening_results: Optional[Dict[str, Any]] = None,
        qm_results: Optional[Dict[str, float]] = None,
        md_results: Optional[Dict[str, float]] = None,
    ) -> float:
        """
        Compute reward for a molecule.
        
        Args:
            smiles: SMILES string
            ml_predictions: ML property predictions
            screening_results: Screening results (similarity, etc.)
            qm_results: QM calculation results
            md_results: MD simulation results
        
        Returns:
            Reward score (higher is better)
        """
        reward = 0.0
        weights = self.config["weights"]
        
        # ML predictions
        if ml_predictions:
            for prop, value in ml_predictions.items():
                if prop in weights:
                    # Normalize if needed
                    normalized = self._normalize_property(prop, value) if self.config["normalization"] else value
                    reward += weights[prop] * normalized
        
        # Screening results
        if screening_results:
            if "similarity" in screening_results and "similarity" in weights:
                sim = screening_results["similarity"]
                reward += weights["similarity"] * sim
        
        # QM results
        if qm_results:
            if "energy" in qm_results and "energy" in weights:
                # Energy is negative (lower is better)
                energy = qm_results["energy"]
                # Normalize: assume energy range [-1000, 0] -> [0, 1]
                normalized_energy = max(0, min(1, (energy + 1000) / 1000))
                reward += weights["energy"] * (1 - normalized_energy)  # Invert
        
        # MD results
        if md_results:
            # Add MD-based rewards if needed
            pass
        
        # Apply constraints
        reward = self._apply_constraints(smiles, reward)
        
        return reward
    
    def _normalize_property(self, prop: str, value: float) -> float:
        """
        Normalize property value to [0, 1] range.
        
        Args:
            prop: Property name
            value: Raw value
        
        Returns:
            Normalized value
        """
        # Mock normalization ranges
        ranges = {
            "logP": (-5, 5),  # Typical logP range
            "solubility": (0, 10),  # -logS range
            "toxicity": (0, 1),  # Binary or probability
        }
        
        if prop in ranges:
            min_val, max_val = ranges[prop]
            normalized = (value - min_val) / (max_val - min_val) if max_val > min_val else 0.5
            return max(0, min(1, normalized))
        
        # Default: assume already in [0, 1]
        return max(0, min(1, value))
    
    def _apply_constraints(self, smiles: str, reward: float) -> float:
        """
        Apply constraint penalties to reward.
        
        Args:
            smiles: SMILES string
            reward: Current reward
        
        Returns:
            Adjusted reward
        """
        constraints = self.config.get("constraints", {})
        
        # Mock: apply simple constraints
        # In real implementation, would check:
        # - Molecular weight limits
        # - Drug-likeness rules (Lipinski, etc.)
        # - Forbidden substructures
        # - Aromatic ring counts
        
        # Example: penalize very long molecules
        if len(smiles) > 100:
            reward *= 0.5
        
        return reward
    
    def batch_compute(
        self,
        smiles_list: List[str],
        predictions_list: Optional[List[Dict[str, float]]] = None,
        screening_list: Optional[List[Dict[str, Any]]] = None,
        qm_list: Optional[List[Dict[str, float]]] = None,
        md_list: Optional[List[Dict[str, float]]] = None,
    ) -> List[float]:
        """
        Compute rewards for a batch of molecules.
        
        Args:
            smiles_list: List of SMILES
            predictions_list: List of ML predictions
            screening_list: List of screening results
            qm_list: List of QM results
            md_list: List of MD results
        
        Returns:
            List of reward scores
        """
        rewards = []
        for i, smiles in enumerate(smiles_list):
            ml_pred = predictions_list[i] if predictions_list and i < len(predictions_list) else None
            screen = screening_list[i] if screening_list and i < len(screening_list) else None
            qm = qm_list[i] if qm_list and i < len(qm_list) else None
            md = md_list[i] if md_list and i < len(md_list) else None
            
            reward = self.compute(smiles, ml_pred, screen, qm, md)
            rewards.append(reward)
        
        return rewards
    
    def update_weights(self, new_weights: Dict[str, float]):
        """
        Update property weights.
        
        Args:
            new_weights: New weight dictionary
        """
        self.config["weights"].update(new_weights)
        logger.info(f"Updated weights: {new_weights}")
    
    def get_config(self) -> Dict[str, Any]:
        """Get reward function configuration."""
        return self.config.copy()

