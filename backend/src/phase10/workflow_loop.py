"""
Workflow Loop - Batch generation → evaluation → policy update

Main loop for RL-based molecule optimization:
1. Generate batch of molecules
2. Evaluate molecules (ML, screening, QM/MD)
3. Compute rewards
4. Update policy
5. Track top candidates
"""

from typing import List, Dict, Optional, Any
import logging
import asyncio

logger = logging.getLogger(__name__)


class WorkflowLoop:
    """
    Main workflow loop for RL-based molecule optimization.
    
    Iteratively:
    1. Generates molecules (RL or generative agent)
    2. Evaluates molecules (ML, screening, QM/MD)
    3. Computes rewards
    4. Updates policy
    5. Tracks top candidates
    """
    
    def __init__(
        self,
        rl_agent: Any,  # RLAgent
        generative_agent: Any,  # GenerativeAgent
        evaluator: Any,  # Evaluator
        dataset_utils: Any,  # DatasetUtils
        config: Optional[Dict[str, Any]] = None,
    ):
        """
        Initialize workflow loop.
        
        Args:
            rl_agent: RLAgent instance
            generative_agent: GenerativeAgent instance
            evaluator: Evaluator instance
            dataset_utils: DatasetUtils instance
            config: Configuration dict with:
                - batch_size: Molecules per iteration
                - max_iterations: Maximum iterations
                - use_generative: Whether to use generative agent
                - top_k: Number of top candidates to track
        """
        self.rl_agent = rl_agent
        self.generative_agent = generative_agent
        self.evaluator = evaluator
        self.dataset_utils = dataset_utils
        self.config = config or {
            "batch_size": 32,
            "max_iterations": 100,
            "use_generative": False,
            "top_k": 10,
        }
        self.iteration_logs: List[Dict[str, Any]] = []
        logger.info("Workflow Loop initialized")
    
    async def run_iteration(
        self,
        iteration: int,
        seed_smiles: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        """
        Run a single iteration of the workflow loop.
        
        Args:
            iteration: Iteration number
            seed_smiles: Optional seed SMILES for generation
        
        Returns:
            Iteration results
        """
        logger.info(f"Starting iteration {iteration}")
        
        # Step 1: Generate molecules
        if self.config.get("use_generative"):
            generated = self.generative_agent.generate(
                n=self.config["batch_size"],
                seed_smiles=seed_smiles,
            )
            generation_method = "generative"
        else:
            generated = self.rl_agent.generate_batch(
                n=self.config["batch_size"],
                seed_smiles=seed_smiles,
            )
            generation_method = "rl"
        
        logger.info(f"Generated {len(generated)} molecules ({generation_method})")
        
        # Step 2: Evaluate molecules
        evaluation_results = await self.evaluator.evaluate_batch(
            generated,
            compute_ml=True,
            compute_screening=True,
            compute_qm=False,  # Expensive, optional
            compute_md=False,  # Expensive, optional
        )
        
        # Step 3: Extract rewards and properties
        rewards = [r.get("reward", 0.0) for r in evaluation_results]
        molecules_with_rewards = [
            (r["smiles"], r["reward"], r.get("ml_predictions", {}))
            for r in evaluation_results
        ]
        
        # Step 4: Update policy (if using RL)
        policy_update = None
        if generation_method == "rl":
            policy_update = self.rl_agent.update_policy(rewards, generated)
            logger.info(f"Policy updated: {policy_update}")
        
        # Step 5: Store in dataset
        for smiles, reward, properties in molecules_with_rewards:
            self.dataset_utils.add_record(
                smiles=smiles,
                reward=reward,
                properties=properties,
                generation_iteration=iteration,
                generation_method=generation_method,
            )
        
        # Step 6: Get top candidates
        top_k = self.config.get("top_k", 10)
        top_candidates = self.dataset_utils.get_top_candidates(n=top_k)
        
        # Step 7: Log iteration
        iteration_log = {
            "iteration": iteration,
            "generation_method": generation_method,
            "generated_count": len(generated),
            "evaluated_count": len(evaluation_results),
            "avg_reward": sum(rewards) / len(rewards) if rewards else 0.0,
            "max_reward": max(rewards) if rewards else 0.0,
            "min_reward": min(rewards) if rewards else 0.0,
            "policy_update": policy_update,
            "top_reward": top_candidates[0].reward if top_candidates else 0.0,
        }
        self.iteration_logs.append(iteration_log)
        
        logger.info(
            f"Iteration {iteration} complete: "
            f"avg_reward={iteration_log['avg_reward']:.4f}, "
            f"max_reward={iteration_log['max_reward']:.4f}"
        )
        
        return iteration_log
    
    async def run(
        self,
        max_iterations: Optional[int] = None,
        seed_smiles: Optional[List[str]] = None,
        early_stop_threshold: Optional[float] = None,
    ) -> Dict[str, Any]:
        """
        Run the full workflow loop.
        
        Args:
            max_iterations: Maximum iterations (overrides config)
            seed_smiles: Optional seed SMILES for first iteration
            early_stop_threshold: Stop if top reward exceeds this
        
        Returns:
            Final results with top candidates and logs
        """
        max_iter = max_iterations or self.config["max_iterations"]
        logger.info(f"Starting workflow loop: {max_iter} iterations")
        
        current_seeds = seed_smiles
        
        for iteration in range(1, max_iter + 1):
            iteration_result = await self.run_iteration(iteration, current_seeds)
            
            # Early stopping
            if early_stop_threshold and iteration_result["max_reward"] >= early_stop_threshold:
                logger.info(f"Early stopping at iteration {iteration} (reward >= {early_stop_threshold})")
                break
            
            # Update seeds for next iteration (use top candidates)
            if iteration_result["max_reward"] > 0:
                top_candidates = self.dataset_utils.get_top_candidates(n=5)
                current_seeds = [c.smiles for c in top_candidates]
        
        # Final summary
        top_candidates = self.dataset_utils.get_top_candidates(n=self.config["top_k"])
        stats = self.dataset_utils.get_statistics()
        
        return {
            "iterations_completed": len(self.iteration_logs),
            "top_candidates": [
                {
                    "smiles": c.smiles,
                    "reward": c.reward,
                    "properties": c.properties,
                    "iteration": c.generation_iteration,
                }
                for c in top_candidates
            ],
            "statistics": stats,
            "iteration_logs": self.iteration_logs,
        }
    
    def get_iteration_logs(self) -> List[Dict[str, Any]]:
        """Get all iteration logs."""
        return self.iteration_logs.copy()
    
    def get_top_candidates(self, n: Optional[int] = None) -> List[Any]:
        """Get top candidates."""
        k = n or self.config["top_k"]
        return self.dataset_utils.get_top_candidates(n=k)

