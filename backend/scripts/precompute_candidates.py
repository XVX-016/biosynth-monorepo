"""
Step 6: Precompute Candidates for Dashboard

Precomputes top molecules for frontend display.
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
)

# Try to import orchestrator
try:
    from src.orchestrator import Orchestrator
    ORCHESTRATOR_AVAILABLE = True
except ImportError:
    ORCHESTRATOR_AVAILABLE = False
    Orchestrator = None

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


async def precompute_candidates(
    batch_size: int = 32,
    iterations: int = 10,
    top_k: int = 20,
):
    """
    Precompute top candidates for frontend dashboard.
    
    Args:
        batch_size: Molecules per iteration
        iterations: Number of iterations
        top_k: Number of top candidates to save
    """
    logger.info("=" * 60)
    logger.info("Precomputing Candidates for Dashboard")
    logger.info("=" * 60)
    
    # Initialize components
    rl_agent = RLAgent()
    generative_agent = GenerativeAgent()
    reward_function = RewardFunction()
    dataset_utils = DatasetUtils(storage_path="data/phase10")
    
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
        generative_agent=generative_agent,
        evaluator=evaluator,
        dataset_utils=dataset_utils,
        config={
            "batch_size": batch_size,
            "max_iterations": iterations,
            "use_generative": False,
            "top_k": top_k,
        },
    )
    
    logger.info(f"\nRunning precomputation:")
    logger.info(f"  Batch size: {batch_size}")
    logger.info(f"  Iterations: {iterations}")
    logger.info(f"  Top K: {top_k}")
    
    # Run workflow
    results = await workflow_loop.run(max_iterations=iterations)
    
    # Get top candidates
    top_candidates = workflow_loop.get_top_candidates(n=top_k)
    
    # Prepare data for frontend
    frontend_data = {
        "timestamp": datetime.now().isoformat(),
        "top_candidates": [
            {
                "smiles": c.smiles,
                "reward": c.reward,
                "properties": c.properties,
                "iteration": c.generation_iteration,
                "method": c.generation_method,
            }
            for c in top_candidates
        ],
        "statistics": results.get("statistics", {}),
        "iteration_logs": results.get("iteration_logs", []),
    }
    
    # Save to frontend cache
    cache_dir = Path("data/frontend_cache")
    cache_dir.mkdir(parents=True, exist_ok=True)
    
    cache_file = cache_dir / "top_candidates.json"
    with open(cache_file, "w") as f:
        json.dump(frontend_data, f, indent=2, default=str)
    
    logger.info(f"\nPrecomputation complete!")
    logger.info(f"  Top candidates: {len(top_candidates)}")
    logger.info(f"  Cache file: {cache_file}")
    
    return frontend_data


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--iterations", type=int, default=10)
    parser.add_argument("--top-k", type=int, default=20)
    
    args = parser.parse_args()
    
    asyncio.run(precompute_candidates(
        batch_size=args.batch_size,
        iterations=args.iterations,
        top_k=args.top_k,
    ))

