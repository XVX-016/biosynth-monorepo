"""
Evaluator - Orchestrates agents to evaluate and score generated molecules

Coordinates ML predictions, screening, QM/MD calculations to compute
rewards for generated molecules.
"""

from typing import List, Dict, Optional, Any
import logging
import asyncio

logger = logging.getLogger(__name__)


class Evaluator:
    """
    Evaluator for generated molecules.
    
    Orchestrates:
    - ML property predictions (Phase 5)
    - Screening operations (Phase 7)
    - QM/MD calculations (Phase 8)
    - Reward computation (Phase 10)
    """
    
    def __init__(
        self,
        reward_function: Optional[Any] = None,  # RewardFunction
        orchestrator: Optional[Any] = None,  # Orchestrator from Phase 9
    ):
        """
        Initialize evaluator.
        
        Args:
            reward_function: RewardFunction instance
            orchestrator: Phase 9 Orchestrator instance
        """
        self.reward_function = reward_function
        self.orchestrator = orchestrator
        logger.info("Evaluator initialized")
    
    async def evaluate_molecule(
        self,
        smiles: str,
        compute_ml: bool = True,
        compute_screening: bool = True,
        compute_qm: bool = False,
        compute_md: bool = False,
    ) -> Dict[str, Any]:
        """
        Evaluate a single molecule.
        
        Args:
            smiles: SMILES string
            compute_ml: Whether to compute ML predictions
            compute_screening: Whether to run screening
            compute_qm: Whether to run QM calculations
            compute_md: Whether to run MD simulations
        
        Returns:
            Evaluation results with reward
        """
        logger.info(f"Evaluating molecule: {smiles}")
        
        results = {
            "smiles": smiles,
            "ml_predictions": {},
            "screening_results": {},
            "qm_results": {},
            "md_results": {},
            "reward": 0.0,
        }
        
        # ML predictions (Phase 5)
        if compute_ml and self.orchestrator:
            try:
                ml_task = {
                    "task_type": "predict",
                    "input_data": {"smiles": smiles, "properties": ["logP", "solubility", "toxicity"]},
                }
                ml_result = await self.orchestrator.execute_task(ml_task)
                if ml_result.success and "output_data" in ml_result.output_data:
                    results["ml_predictions"] = ml_result.output_data.get("predictions", {})
            except Exception as e:
                logger.warning(f"ML prediction failed: {e}")
        
        # Screening (Phase 7)
        if compute_screening and self.orchestrator:
            try:
                screen_task = {
                    "task_type": "screen",
                    "input_data": {
                        "task_subtype": "similarity",
                        "query_smiles": smiles,
                        "k": 5,
                    },
                }
                screen_result = await self.orchestrator.execute_task(screen_task)
                if screen_result.success and "output_data" in screen_result.output_data:
                    results["screening_results"] = screen_result.output_data.get("results", {})
            except Exception as e:
                logger.warning(f"Screening failed: {e}")
        
        # QM calculations (Phase 8)
        if compute_qm and self.orchestrator:
            try:
                qm_task = {
                    "task_type": "qm",
                    "input_data": {"smiles": smiles, "method": "energy"},
                }
                qm_result = await self.orchestrator.execute_task(qm_task)
                if qm_result.success and "output_data" in qm_result.output_data:
                    results["qm_results"] = qm_result.output_data.get("results", {})
            except Exception as e:
                logger.warning(f"QM calculation failed: {e}")
        
        # MD simulations (Phase 8)
        if compute_md and self.orchestrator:
            try:
                md_task = {
                    "task_type": "md",
                    "input_data": {"smiles": smiles, "steps": 100},
                }
                md_result = await self.orchestrator.execute_task(md_task)
                if md_result.success and "output_data" in md_result.output_data:
                    results["md_results"] = md_result.output_data.get("results", {})
            except Exception as e:
                logger.warning(f"MD simulation failed: {e}")
        
        # Compute reward
        if self.reward_function:
            results["reward"] = self.reward_function.compute(
                smiles,
                results["ml_predictions"],
                results["screening_results"],
                results["qm_results"],
                results["md_results"],
            )
        
        return results
    
    async def evaluate_batch(
        self,
        smiles_list: List[str],
        compute_ml: bool = True,
        compute_screening: bool = True,
        compute_qm: bool = False,
        compute_md: bool = False,
    ) -> List[Dict[str, Any]]:
        """
        Evaluate a batch of molecules.
        
        Args:
            smiles_list: List of SMILES
            compute_ml: Whether to compute ML predictions
            compute_screening: Whether to run screening
            compute_qm: Whether to run QM calculations
            compute_md: Whether to run MD simulations
        
        Returns:
            List of evaluation results
        """
        logger.info(f"Evaluating batch of {len(smiles_list)} molecules")
        
        # Evaluate in parallel (with limit to avoid overwhelming)
        tasks = [
            self.evaluate_molecule(
                smiles,
                compute_ml=compute_ml,
                compute_screening=compute_screening,
                compute_qm=compute_qm,
                compute_md=compute_md,
            )
            for smiles in smiles_list
        ]
        
        results = await asyncio.gather(*tasks, return_exceptions=True)
        
        # Filter out exceptions
        valid_results = []
        for r in results:
            if isinstance(r, Exception):
                logger.error(f"Evaluation error: {r}")
            else:
                valid_results.append(r)
        
        return valid_results
    
    def get_evaluation_summary(self, results: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Get summary statistics from evaluation results.
        
        Args:
            results: List of evaluation results
        
        Returns:
            Summary statistics
        """
        if not results:
            return {"count": 0}
        
        rewards = [r.get("reward", 0.0) for r in results]
        
        return {
            "count": len(results),
            "avg_reward": sum(rewards) / len(rewards) if rewards else 0.0,
            "max_reward": max(rewards) if rewards else 0.0,
            "min_reward": min(rewards) if rewards else 0.0,
            "success_rate": len([r for r in results if r.get("reward", 0) > 0]) / len(results),
        }

