"""
Phase D: Orchestrator Stabilization

Goal: Ensure Phase 9 agent routing handles failures and edge cases.

Tasks:
1. Introduce failures (invalid SMILES, failed QM/MD)
2. Verify orchestrator retries tasks and logs failures
3. Ensure dependent tasks continue correctly
"""

import sys
from pathlib import Path
import logging
import asyncio
from typing import Dict, Any, List

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

try:
    from src.orchestrator import Orchestrator, Task, TaskResult
    ORCHESTRATOR_AVAILABLE = True
except ImportError:
    ORCHESTRATOR_AVAILABLE = False
    logging.warning("Orchestrator not available")

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class OrchestratorStabilizationTester:
    """Tests orchestrator error handling and stability."""
    
    def __init__(self):
        self.results = {
            "timestamp": None,
            "failure_tests": [],
            "retry_tests": [],
            "dependent_task_tests": [],
            "errors": [],
        }
    
    async def test_failure_handling(self) -> Dict[str, Any]:
        """Test orchestrator handles failures gracefully."""
        if not ORCHESTRATOR_AVAILABLE:
            logger.warning("Orchestrator not available. Skipping tests.")
            return {}
        
        logger.info("=" * 60)
        logger.info("Phase D: Orchestrator Stabilization")
        logger.info("=" * 60)
        
        orchestrator = Orchestrator()
        
        # Test 1: Invalid SMILES
        logger.info("\nTest 1: Invalid SMILES handling")
        try:
            task = Task(
                task_type="predict",
                input_data={"smiles": "INVALID_SMILES_!!!", "properties": ["logP"]},
            )
            result = await orchestrator.execute_task(task)
            
            if result.success:
                logger.warning("  ⚠ Task succeeded with invalid SMILES (unexpected)")
            else:
                logger.info(f"  ✓ Task correctly failed: {result.error}")
            
            self.results["failure_tests"].append({
                "test": "invalid_smiles",
                "success": not result.success,
                "error": result.error,
            })
        except Exception as e:
            logger.error(f"  ✗ Test failed: {e}")
            self.results["errors"].append({"test": "invalid_smiles", "error": str(e)})
        
        # Test 2: Missing required fields
        logger.info("\nTest 2: Missing required fields")
        try:
            task = Task(
                task_type="predict",
                input_data={},  # Missing smiles
            )
            result = await orchestrator.execute_task(task)
            
            if result.success:
                logger.warning("  ⚠ Task succeeded with missing fields (unexpected)")
            else:
                logger.info(f"  ✓ Task correctly failed: {result.error}")
            
            self.results["failure_tests"].append({
                "test": "missing_fields",
                "success": not result.success,
                "error": result.error,
            })
        except Exception as e:
            logger.error(f"  ✗ Test failed: {e}")
            self.results["errors"].append({"test": "missing_fields", "error": str(e)})
        
        # Test 3: Unknown task type
        logger.info("\nTest 3: Unknown task type")
        try:
            task = Task(
                task_type="unknown_task_type",
                input_data={"smiles": "CCO"},
            )
            result = await orchestrator.execute_task(task)
            
            if result.success:
                logger.warning("  ⚠ Task succeeded with unknown type (unexpected)")
            else:
                logger.info(f"  ✓ Task correctly failed: {result.error}")
            
            self.results["failure_tests"].append({
                "test": "unknown_task_type",
                "success": not result.success,
                "error": result.error,
            })
        except Exception as e:
            logger.error(f"  ✗ Test failed: {e}")
            self.results["errors"].append({"test": "unknown_task_type", "error": str(e)})
        
        return self.results
    
    async def test_dependent_tasks(self) -> Dict[str, Any]:
        """Test dependent task execution."""
        if not ORCHESTRATOR_AVAILABLE:
            return {}
        
        logger.info("\n" + "=" * 60)
        logger.info("Dependent Task Tests")
        logger.info("=" * 60)
        
        orchestrator = Orchestrator()
        
        # Test: Workflow with one failing task
        logger.info("\nTest: Workflow with failing task")
        try:
            workflow = [
                {
                    "task_type": "predict",
                    "input_data": {"smiles": "CCO", "properties": ["logP"]},
                },
                {
                    "task_type": "predict",
                    "input_data": {"smiles": "INVALID!!!", "properties": ["logP"]},  # This will fail
                },
                {
                    "task_type": "predict",
                    "input_data": {"smiles": "CCCO", "properties": ["logP"]},
                },
            ]
            
            results = await orchestrator.execute_workflow(workflow)
            
            success_count = sum(1 for r in results if r.success)
            logger.info(f"  Workflow completed: {success_count}/{len(results)} tasks successful")
            
            self.results["dependent_task_tests"].append({
                "test": "workflow_with_failure",
                "total_tasks": len(results),
                "successful_tasks": success_count,
            })
        except Exception as e:
            logger.error(f"  ✗ Test failed: {e}")
            self.results["errors"].append({"test": "workflow_with_failure", "error": str(e)})
        
        return self.results
    
    def save_results(self, output_file: str = "data/validation/phase_d_results.json"):
        """Save stabilization test results."""
        from datetime import datetime
        self.results["timestamp"] = datetime.now().isoformat()
        
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        import json
        with open(output_path, "w") as f:
            json.dump(self.results, f, indent=2, default=str)
        
        logger.info(f"\nResults saved to: {output_path}")


async def main():
    """Run Phase D stabilization tests."""
    tester = OrchestratorStabilizationTester()
    
    # Test failure handling
    await tester.test_failure_handling()
    
    # Test dependent tasks
    await tester.test_dependent_tasks()
    
    # Save results
    tester.save_results()
    
    logger.info("\n" + "=" * 60)
    logger.info("Stabilization test complete!")
    logger.info("=" * 60)


if __name__ == "__main__":
    asyncio.run(main())

