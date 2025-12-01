"""
Run All Post-Phase 10 Validation Phases

Orchestrates execution of all validation phases (A-G).
"""

import sys
from pathlib import Path
import logging
import asyncio

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

from scripts.phase_a_end_to_end_validation import EndToEndValidator
from scripts.phase_b_rl_sanity_check import RLSanityChecker
from scripts.phase_c_frontend_integration_test import FrontendIntegrationTester
from scripts.phase_d_orchestrator_stabilization import OrchestratorStabilizationTester
from scripts.phase_e_baseline_metrics import BaselineMetricsCollector
from scripts.phase_f_documentation import DocumentationGenerator
from scripts.phase_g_scaling_decision import ScalingDecisionFramework

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


async def run_all_phases():
    """Run all post-phase 10 validation phases."""
    logger.info("=" * 80)
    logger.info("Post-Phase 10 Validation - Complete Suite")
    logger.info("=" * 80)
    
    results = {}
    
    # Phase A: End-to-End Validation
    logger.info("\n" + "=" * 80)
    logger.info("PHASE A: End-to-End Validation")
    logger.info("=" * 80)
    try:
        validator = EndToEndValidator()
        test_molecules = ["CCO", "CCCO", "CC(C)O", "c1ccccc1", "CCc1ccccc1"]
        results["phase_a"] = await validator.validate_pipeline(test_molecules)
        validator.save_results()
    except Exception as e:
        logger.error(f"Phase A failed: {e}")
        results["phase_a"] = {"error": str(e)}
    
    # Phase B: RL Sanity Check
    logger.info("\n" + "=" * 80)
    logger.info("PHASE B: RL Loop Sanity Check")
    logger.info("=" * 80)
    try:
        checker = RLSanityChecker()
        results["phase_b"] = await checker.run_sanity_check(
            batch_size=5,
            iterations=5,
            seed_smiles=["CCO", "CCCO"],
        )
        checker.save_results()
    except Exception as e:
        logger.error(f"Phase B failed: {e}")
        results["phase_b"] = {"error": str(e)}
    
    # Phase C: Frontend Integration (API tests only)
    logger.info("\n" + "=" * 80)
    logger.info("PHASE C: Backend → Frontend Integration")
    logger.info("=" * 80)
    try:
        tester = FrontendIntegrationTester()
        tester.test_api_endpoints()
        tester.test_data_flow()
        results["phase_c"] = tester.results
        tester.save_results()
    except Exception as e:
        logger.error(f"Phase C failed: {e}")
        results["phase_c"] = {"error": str(e)}
    
    # Phase D: Orchestrator Stabilization
    logger.info("\n" + "=" * 80)
    logger.info("PHASE D: Orchestrator Stabilization")
    logger.info("=" * 80)
    try:
        stabilizer = OrchestratorStabilizationTester()
        await stabilizer.test_failure_handling()
        await stabilizer.test_dependent_tasks()
        results["phase_d"] = stabilizer.results
        stabilizer.save_results()
    except Exception as e:
        logger.error(f"Phase D failed: {e}")
        results["phase_d"] = {"error": str(e)}
    
    # Phase E: Baseline Metrics
    logger.info("\n" + "=" * 80)
    logger.info("PHASE E: Baseline Metrics Collection")
    logger.info("=" * 80)
    try:
        collector = BaselineMetricsCollector()
        results["phase_e"] = await collector.collect_all_metrics()
        collector.save_metrics()
    except Exception as e:
        logger.error(f"Phase E failed: {e}")
        results["phase_e"] = {"error": str(e)}
    
    # Phase F: Documentation
    logger.info("\n" + "=" * 80)
    logger.info("PHASE F: Documentation & Reproducibility")
    logger.info("=" * 80)
    try:
        doc_generator = DocumentationGenerator()
        results["phase_f"] = doc_generator.generate_documentation()
        doc_generator.save_documentation()
    except Exception as e:
        logger.error(f"Phase F failed: {e}")
        results["phase_f"] = {"error": str(e)}
    
    # Phase G: Scaling Decision
    logger.info("\n" + "=" * 80)
    logger.info("PHASE G: Scaling Decision Framework")
    logger.info("=" * 80)
    try:
        framework = ScalingDecisionFramework()
        results["phase_g"] = framework.analyze_scaling()
        framework.save_analysis()
    except Exception as e:
        logger.error(f"Phase G failed: {e}")
        results["phase_g"] = {"error": str(e)}
    
    # Summary
    logger.info("\n" + "=" * 80)
    logger.info("VALIDATION SUITE COMPLETE")
    logger.info("=" * 80)
    
    success_count = sum(1 for phase, result in results.items() if "error" not in result)
    total_phases = len(results)
    
    logger.info(f"\nPhases completed: {success_count}/{total_phases}")
    logger.info("\nResults saved to:")
    logger.info("  - data/validation/phase_a_results.json")
    logger.info("  - data/validation/phase_b_results.json")
    logger.info("  - data/validation/phase_c_results.json")
    logger.info("  - data/validation/phase_d_results.json")
    logger.info("  - data/metrics/baseline_metrics.json")
    logger.info("  - data/docs/reproducibility_docs.json")
    logger.info("  - data/docs/scaling_analysis.json")
    
    return results


if __name__ == "__main__":
    asyncio.run(run_all_phases())

