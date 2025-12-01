"""
Step 1: Run Small RL Loop Test

Executes Phase 10 RL loop with small test batch.
"""

import sys
from pathlib import Path
import logging
import asyncio
import json
from datetime import datetime

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
    Phase10Orchestrator,
)

# Try to import orchestrator for integration
try:
    from src.orchestrator import Orchestrator
    ORCHESTRATOR_AVAILABLE = True
except ImportError:
    ORCHESTRATOR_AVAILABLE = False
    Orchestrator = None

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


async def run_rl_loop_test():
    """
    Run small RL loop test.
    
    Uses 5-10 molecules per iteration, 5 iterations.
    """
    logger.info("=" * 60)
    logger.info("Phase 10 RL Loop Test")
    logger.info("=" * 60)
    
    # Initialize components
    rl_agent = RLAgent(config={
        "learning_rate": 0.001,
        "batch_size": 5,  # Small batch for testing
        "exploration_rate": 0.1,
    })
    
    generative_agent = GenerativeAgent()
    reward_function = RewardFunction()
    dataset_utils = DatasetUtils(storage_path="data/phase10")
    
    # Initialize orchestrator if available
    orchestrator = None
    if ORCHESTRATOR_AVAILABLE:
        try:
            orchestrator = Orchestrator()
            logger.info("Orchestrator initialized")
        except Exception as e:
            logger.warning(f"Could not initialize orchestrator: {e}")
    
    evaluator = Evaluator(
        reward_function=reward_function,
        orchestrator=orchestrator,
    )
    
    workflow_loop = WorkflowLoop(
        rl_agent=rl_agent,
        generative_agent=generative_agent,
        evaluator=evaluator,
        dataset_utils=dataset_utils,
        config={
            "batch_size": 5,
            "max_iterations": 5,
            "use_generative": False,
            "top_k": 10,
        },
    )
    
    # Seed molecules (optional)
    seed_smiles = ["CCO", "CCCO", "CC(C)O"]
    
    logger.info(f"\nStarting RL loop:")
    logger.info(f"  Batch size: 5")
    logger.info(f"  Iterations: 5")
    logger.info(f"  Seed molecules: {len(seed_smiles)}")
    
    # Run workflow
    try:
        results = await workflow_loop.run(
            max_iterations=5,
            seed_smiles=seed_smiles,
        )
        
        logger.info("\n" + "=" * 60)
        logger.info("RL Loop Complete!")
        logger.info("=" * 60)
        logger.info(f"Iterations completed: {results['iterations_completed']}")
        logger.info(f"Total molecules generated: {results['statistics'].get('total', 0)}")
        logger.info(f"Top reward: {results['top_candidates'][0]['reward']:.4f}" if results['top_candidates'] else "N/A")
        
        # Display top candidates
        logger.info("\nTop 5 Candidates:")
        for i, candidate in enumerate(results['top_candidates'][:5], 1):
            logger.info(f"  {i}. {candidate['smiles']} - Reward: {candidate['reward']:.4f}")
        
        # Save results
        output_dir = Path("data/phase10/rl_logs")
        output_dir.mkdir(parents=True, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = output_dir / f"rl_test_{timestamp}.json"
        
        with open(output_file, "w") as f:
            json.dump(results, f, indent=2, default=str)
        
        logger.info(f"\nResults saved to: {output_file}")
        
        return results
        
    except Exception as e:
        logger.error(f"RL loop failed: {e}")
        import traceback
        traceback.print_exc()
        return None


if __name__ == "__main__":
    asyncio.run(run_rl_loop_test())

