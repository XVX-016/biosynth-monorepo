"""
Phase B: RL Loop Sanity Check

Goal: Validate reward function, generative output, and loop convergence.

Tasks:
1. Run a micro RL batch (5-10 molecules)
2. Record rewards at each iteration
3. Detect invalid SMILES or failed conformers
4. Confirm RL candidates improve metrics over baseline
"""

import sys
from pathlib import Path
import logging
import asyncio
import json
from datetime import datetime
from typing import List, Dict, Any

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

from src.phase10 import (
    RLAgent,
    GenerativeAgent,
    RewardFunction,
    Evaluator,
    WorkflowLoop,
    DatasetUtils,
)

try:
    from src.orchestrator import Orchestrator
    ORCHESTRATOR_AVAILABLE = True
except ImportError:
    ORCHESTRATOR_AVAILABLE = False
    Orchestrator = None

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class RLSanityChecker:
    """Validates RL loop sanity and convergence."""
    
    def __init__(self):
        self.results = {
            "timestamp": datetime.now().isoformat(),
            "iterations": [],
            "reward_progression": [],
            "invalid_smiles": [],
            "failed_conformers": [],
            "improvement_detected": False,
        }
    
    async def run_sanity_check(
        self,
        batch_size: int = 5,
        iterations: int = 5,
        seed_smiles: List[str] = None,
    ) -> Dict[str, Any]:
        """
        Run RL loop sanity check.
        
        Args:
            batch_size: Molecules per iteration
            iterations: Number of iterations
            seed_smiles: Optional seed molecules
        """
        logger.info("=" * 60)
        logger.info("Phase B: RL Loop Sanity Check")
        logger.info("=" * 60)
        
        # Initialize components
        rl_agent = RLAgent(config={"batch_size": batch_size})
        gen_agent = GenerativeAgent()
        reward_function = RewardFunction()
        dataset_utils = DatasetUtils()
        
        orchestrator = None
        if ORCHESTRATOR_AVAILABLE:
            try:
                orchestrator = Orchestrator()
            except Exception as e:
                logger.warning(f"Could not initialize orchestrator: {e}")
        
        evaluator = Evaluator(
            reward_function=reward_function,
            orchestrator=orchestrator,
        )
        
        workflow_loop = WorkflowLoop(
            rl_agent=rl_agent,
            generative_agent=gen_agent,
            evaluator=evaluator,
            dataset_utils=dataset_utils,
            config={
                "batch_size": batch_size,
                "max_iterations": iterations,
                "use_generative": False,
                "top_k": 10,
            },
        )
        
        # Baseline: evaluate seed molecules
        baseline_rewards = []
        if seed_smiles:
            logger.info(f"\nEvaluating baseline molecules: {seed_smiles}")
            for smiles in seed_smiles:
                try:
                    result = await evaluator.evaluate_molecule(smiles)
                    baseline_rewards.append(result["reward"])
                except Exception as e:
                    logger.warning(f"Baseline evaluation failed for {smiles}: {e}")
        
        baseline_avg = sum(baseline_rewards) / len(baseline_rewards) if baseline_rewards else 0.0
        logger.info(f"Baseline average reward: {baseline_avg:.4f}")
        
        # Run RL loop
        logger.info(f"\nRunning RL loop: {iterations} iterations, {batch_size} molecules each")
        
        try:
            results = await workflow_loop.run(
                max_iterations=iterations,
                seed_smiles=seed_smiles,
            )
            
            # Analyze results
            self._analyze_results(results, baseline_avg)
            
            self.results["rl_results"] = results
            self.results["baseline_avg"] = baseline_avg
            
        except Exception as e:
            logger.error(f"RL loop failed: {e}")
            self.results["error"] = str(e)
        
        return self.results
    
    def _analyze_results(self, results: Dict[str, Any], baseline_avg: float):
        """Analyze RL loop results."""
        iteration_logs = results.get("iteration_logs", [])
        top_candidates = results.get("top_candidates", [])
        
        logger.info("\n" + "=" * 60)
        logger.info("RL Loop Analysis")
        logger.info("=" * 60)
        
        # Reward progression
        rewards_by_iteration = []
        for log in iteration_logs:
            iter_num = log.get("iteration", 0)
            avg_reward = log.get("avg_reward", 0.0)
            max_reward = log.get("max_reward", 0.0)
            
            rewards_by_iteration.append({
                "iteration": iter_num,
                "avg_reward": avg_reward,
                "max_reward": max_reward,
            })
            
            logger.info(f"Iteration {iter_num}: avg={avg_reward:.4f}, max={max_reward:.4f}")
        
        self.results["reward_progression"] = rewards_by_iteration
        
        # Check for improvement
        if rewards_by_iteration:
            final_max = rewards_by_iteration[-1]["max_reward"]
            initial_max = rewards_by_iteration[0]["max_reward"] if rewards_by_iteration else 0.0
            
            improvement = final_max > baseline_avg and final_max > initial_max
            self.results["improvement_detected"] = improvement
            
            if improvement:
                logger.info(f"\n✓ Improvement detected!")
                logger.info(f"  Baseline: {baseline_avg:.4f}")
                logger.info(f"  Initial max: {initial_max:.4f}")
                logger.info(f"  Final max: {final_max:.4f}")
            else:
                logger.warning(f"\n⚠ No clear improvement detected")
                logger.warning(f"  Baseline: {baseline_avg:.4f}")
                logger.warning(f"  Final max: {final_max:.4f}")
        
        # Top candidates
        if top_candidates:
            logger.info(f"\nTop 3 Candidates:")
            for i, candidate in enumerate(top_candidates[:3], 1):
                logger.info(f"  {i}. {candidate['smiles']} - Reward: {candidate['reward']:.4f}")
    
    def save_results(self, output_file: str = "data/validation/phase_b_results.json"):
        """Save sanity check results."""
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, "w") as f:
            json.dump(self.results, f, indent=2, default=str)
        
        logger.info(f"\nResults saved to: {output_path}")


async def main():
    """Run Phase B sanity check."""
    checker = RLSanityChecker()
    
    seed_smiles = ["CCO", "CCCO", "CC(C)O"]
    
    results = await checker.run_sanity_check(
        batch_size=5,
        iterations=5,
        seed_smiles=seed_smiles,
    )
    
    checker.save_results()


if __name__ == "__main__":
    asyncio.run(main())

