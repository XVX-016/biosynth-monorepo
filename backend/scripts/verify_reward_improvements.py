"""
Step 2: Verify Reward Improvements

Checks that RL loop produces improving scores over iterations.
"""

import sys
from pathlib import Path
import logging
import json
from typing import List, Dict, Any

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def verify_reward_improvements(log_file: str = None):
    """
    Verify that RL loop produces improving reward scores.
    
    Args:
        log_file: Path to RL log file (optional, will find latest if not provided)
    """
    logger.info("=" * 60)
    logger.info("Verifying Reward Improvements")
    logger.info("=" * 60)
    
    logs_dir = Path("data/phase10/rl_logs")
    
    if not logs_dir.exists():
        logger.error(f"Logs directory not found: {logs_dir}")
        return False
    
    # Find log file
    if log_file:
        log_path = Path(log_file)
    else:
        # Find latest log file
        log_files = sorted(logs_dir.glob("rl_test_*.json"), reverse=True)
        if not log_files:
            logger.error("No log files found")
            return False
        log_path = log_files[0]
    
    logger.info(f"Analyzing log file: {log_path}")
    
    # Load log
    with open(log_path, "r") as f:
        results = json.load(f)
    
    iteration_logs = results.get("iteration_logs", [])
    
    if not iteration_logs:
        logger.error("No iteration logs found")
        return False
    
    logger.info(f"\nFound {len(iteration_logs)} iterations")
    
    # Extract reward progression
    avg_rewards = []
    max_rewards = []
    iterations = []
    
    for log in iteration_logs:
        iterations.append(log.get("iteration", 0))
        avg_rewards.append(log.get("avg_reward", 0.0))
        max_rewards.append(log.get("max_reward", 0.0))
    
    # Analyze progression
    logger.info("\nReward Progression:")
    logger.info("Iteration | Avg Reward | Max Reward")
    logger.info("-" * 40)
    for i, (iter_num, avg, max_r) in enumerate(zip(iterations, avg_rewards, max_rewards)):
        logger.info(f"    {iter_num:2d}   |   {avg:7.4f}  |  {max_r:7.4f}")
    
    # Check for improvement
    improvements = []
    for i in range(1, len(max_rewards)):
        if max_rewards[i] > max_rewards[i-1]:
            improvements.append(True)
        else:
            improvements.append(False)
    
    improvement_rate = sum(improvements) / len(improvements) if improvements else 0.0
    
    logger.info(f"\nImprovement Analysis:")
    logger.info(f"  Iterations with improvement: {sum(improvements)}/{len(improvements)}")
    logger.info(f"  Improvement rate: {improvement_rate:.2%}")
    logger.info(f"  Final max reward: {max_rewards[-1]:.4f}")
    logger.info(f"  Initial max reward: {max_rewards[0]:.4f}")
    logger.info(f"  Total improvement: {max_rewards[-1] - max_rewards[0]:.4f}")
    
    # Check top candidates
    top_candidates = results.get("top_candidates", [])
    if top_candidates:
        logger.info(f"\nTop 3 Candidates:")
        for i, candidate in enumerate(top_candidates[:3], 1):
            logger.info(f"  {i}. {candidate['smiles']} - Reward: {candidate['reward']:.4f}")
    
    # Summary
    is_improving = improvement_rate > 0.3 or (max_rewards[-1] > max_rewards[0])
    
    if is_improving:
        logger.info("\n✓ RL loop is producing improving rewards")
        return True
    else:
        logger.warning("\n⚠ RL loop may not be improving (check reward function and policy updates)")
        return False


if __name__ == "__main__":
    verify_reward_improvements()

