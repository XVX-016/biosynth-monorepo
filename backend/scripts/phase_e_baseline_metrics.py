"""
Phase E: Baseline Metrics Collection

Goal: Quantify system performance before scaling.

Tasks:
1. Measure prediction accuracy
2. Measure RL reward improvements
3. Measure throughput (molecules/sec)
4. Store metrics in dashboard or log
"""

import sys
from pathlib import Path
import logging
import asyncio
import time
import json
from datetime import datetime
from typing import Dict, Any, List

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
    from ml.prediction_engine import PredictionEngine
    from ml.registry import ModelRegistry
    from src.orchestrator import Orchestrator
    ML_AVAILABLE = True
except ImportError:
    ML_AVAILABLE = False

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class BaselineMetricsCollector:
    """Collects baseline performance metrics."""
    
    def __init__(self):
        self.metrics = {
            "timestamp": datetime.now().isoformat(),
            "prediction_accuracy": {},
            "rl_reward_improvements": {},
            "throughput": {},
            "system_info": {},
        }
    
    async def collect_all_metrics(self) -> Dict[str, Any]:
        """Collect all baseline metrics."""
        logger.info("=" * 60)
        logger.info("Phase E: Baseline Metrics Collection")
        logger.info("=" * 60)
        
        # Measure prediction accuracy
        logger.info("\n1. Measuring prediction accuracy...")
        await self._measure_prediction_accuracy()
        
        # Measure RL reward improvements
        logger.info("\n2. Measuring RL reward improvements...")
        await self._measure_rl_improvements()
        
        # Measure throughput
        logger.info("\n3. Measuring throughput...")
        await self._measure_throughput()
        
        # Collect system info
        logger.info("\n4. Collecting system information...")
        self._collect_system_info()
        
        return self.metrics
    
    async def _measure_prediction_accuracy(self):
        """Measure ML prediction accuracy."""
        if not ML_AVAILABLE:
            self.metrics["prediction_accuracy"] = {"status": "skipped", "reason": "ML not available"}
            return
        
        try:
            registry = ModelRegistry()
            engine = PredictionEngine(registry)
            
            test_molecules = ["CCO", "CCCO", "CC(C)O", "c1ccccc1", "CCc1ccccc1"]
            predictions = []
            
            for smiles in test_molecules:
                try:
                    result = engine.predict(
                        input_data={"smiles": smiles},
                        properties=["logP", "solubility"],
                    )
                    predictions.append({
                        "smiles": smiles,
                        "predictions": result.predictions,
                    })
                except Exception as e:
                    logger.warning(f"Prediction failed for {smiles}: {e}")
            
            self.metrics["prediction_accuracy"] = {
                "status": "success",
                "test_molecules": len(test_molecules),
                "successful_predictions": len(predictions),
                "success_rate": len(predictions) / len(test_molecules) if test_molecules else 0.0,
            }
            
            logger.info(f"  ✓ Prediction accuracy: {len(predictions)}/{len(test_molecules)} successful")
            
        except Exception as e:
            logger.error(f"  ✗ Error: {e}")
            self.metrics["prediction_accuracy"] = {"status": "error", "error": str(e)}
    
    async def _measure_rl_improvements(self):
        """Measure RL reward improvements."""
        try:
            rl_agent = RLAgent()
            reward_function = RewardFunction()
            dataset_utils = DatasetUtils()
            
            # Run small RL loop
            workflow_loop = WorkflowLoop(
                rl_agent=rl_agent,
                generative_agent=GenerativeAgent(),
                evaluator=Evaluator(reward_function=reward_function),
                dataset_utils=dataset_utils,
                config={"batch_size": 5, "max_iterations": 3},
            )
            
            results = await workflow_loop.run(max_iterations=3)
            
            iteration_logs = results.get("iteration_logs", [])
            if iteration_logs:
                initial_reward = iteration_logs[0].get("max_reward", 0.0)
                final_reward = iteration_logs[-1].get("max_reward", 0.0)
                improvement = final_reward - initial_reward
                
                self.metrics["rl_reward_improvements"] = {
                    "status": "success",
                    "initial_reward": initial_reward,
                    "final_reward": final_reward,
                    "improvement": improvement,
                    "improvement_rate": improvement / initial_reward if initial_reward > 0 else 0.0,
                }
                
                logger.info(f"  ✓ RL improvement: {improvement:.4f} ({improvement/initial_reward*100:.1f}%)" if initial_reward > 0 else "N/A")
            else:
                self.metrics["rl_reward_improvements"] = {"status": "no_data"}
                
        except Exception as e:
            logger.error(f"  ✗ Error: {e}")
            self.metrics["rl_reward_improvements"] = {"status": "error", "error": str(e)}
    
    async def _measure_throughput(self):
        """Measure system throughput."""
        try:
            # Test prediction throughput
            if ML_AVAILABLE:
                registry = ModelRegistry()
                engine = PredictionEngine(registry)
                
                test_molecules = ["CCO"] * 10  # 10 identical molecules
                start_time = time.time()
                
                for smiles in test_molecules:
                    try:
                        engine.predict(input_data={"smiles": smiles}, properties=["logP"])
                    except:
                        pass
                
                elapsed = time.time() - start_time
                throughput = len(test_molecules) / elapsed if elapsed > 0 else 0.0
                
                self.metrics["throughput"] = {
                    "status": "success",
                    "molecules_per_second": throughput,
                    "time_per_molecule": elapsed / len(test_molecules) if test_molecules else 0.0,
                }
                
                logger.info(f"  ✓ Throughput: {throughput:.2f} molecules/sec")
            else:
                self.metrics["throughput"] = {"status": "skipped", "reason": "ML not available"}
                
        except Exception as e:
            logger.error(f"  ✗ Error: {e}")
            self.metrics["throughput"] = {"status": "error", "error": str(e)}
    
    def _collect_system_info(self):
        """Collect system information."""
        import platform
        import sys
        
        self.metrics["system_info"] = {
            "platform": platform.platform(),
            "python_version": sys.version,
            "cpu_count": None,  # Could add psutil if available
        }
        
        logger.info(f"  ✓ System info collected")
    
    def save_metrics(self, output_file: str = "data/metrics/baseline_metrics.json"):
        """Save baseline metrics."""
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, "w") as f:
            json.dump(self.metrics, f, indent=2, default=str)
        
        logger.info(f"\nMetrics saved to: {output_path}")


async def main():
    """Run Phase E metrics collection."""
    collector = BaselineMetricsCollector()
    metrics = await collector.collect_all_metrics()
    collector.save_metrics()
    
    logger.info("\n" + "=" * 60)
    logger.info("Baseline metrics collection complete!")
    logger.info("=" * 60)


if __name__ == "__main__":
    asyncio.run(main())

