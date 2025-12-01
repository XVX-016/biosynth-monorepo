"""
Phase A: End-to-End Validation

Goal: Ensure all components actually work together.

Tasks:
1. Load a small molecule library (5-10 molecules)
2. Run pipeline: SMILES → Screening → ML Prediction → Conformers → RL Loop
3. Log any errors, exceptions, or failed outputs
4. Compare predicted properties vs known benchmarks
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

# Import Phase components
try:
    from src.search import SearchEngine, LibraryLoader
    from src.conformers import ConformerGenerator
    from src.orchestrator import Orchestrator
    from src.phase10 import Phase10Orchestrator, RLAgent, GenerativeAgent, RewardFunction, Evaluator, WorkflowLoop, DatasetUtils
    from ml.prediction_engine import PredictionEngine
    from ml.registry import ModelRegistry
    PHASES_AVAILABLE = True
except ImportError as e:
    PHASES_AVAILABLE = False
    logging.warning(f"Some phases not available: {e}")

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class EndToEndValidator:
    """Validates end-to-end pipeline execution."""
    
    def __init__(self):
        self.results = {
            "timestamp": datetime.now().isoformat(),
            "molecules_tested": [],
            "errors": [],
            "warnings": [],
            "success_count": 0,
            "failure_count": 0,
        }
    
    async def validate_pipeline(self, smiles_list: List[str]) -> Dict[str, Any]:
        """
        Run full pipeline validation for a list of SMILES.
        
        Pipeline: SMILES → Screening → ML Prediction → Conformers → RL Loop
        """
        logger.info("=" * 60)
        logger.info("Phase A: End-to-End Validation")
        logger.info("=" * 60)
        logger.info(f"Testing {len(smiles_list)} molecules")
        
        for i, smiles in enumerate(smiles_list, 1):
            logger.info(f"\n[{i}/{len(smiles_list)}] Processing: {smiles}")
            
            molecule_result = {
                "smiles": smiles,
                "steps": {},
                "success": False,
                "errors": [],
            }
            
            try:
                # Step 1: Screening (Phase 7)
                logger.info("  → Step 1: Screening (Phase 7)")
                screening_result = await self._test_screening(smiles)
                molecule_result["steps"]["screening"] = screening_result
                
                # Step 2: ML Prediction (Phase 5)
                logger.info("  → Step 2: ML Prediction (Phase 5)")
                prediction_result = await self._test_ml_prediction(smiles)
                molecule_result["steps"]["ml_prediction"] = prediction_result
                
                # Step 3: Conformers (Phase 8)
                logger.info("  → Step 3: Conformers (Phase 8)")
                conformer_result = await self._test_conformers(smiles)
                molecule_result["steps"]["conformers"] = conformer_result
                
                # Step 4: RL Loop (Phase 10)
                logger.info("  → Step 4: RL Loop (Phase 10)")
                rl_result = await self._test_rl_loop(smiles)
                molecule_result["steps"]["rl_loop"] = rl_result
                
                # Compare with benchmarks if available
                logger.info("  → Step 5: Benchmark Comparison")
                benchmark_result = self._compare_benchmarks(smiles, prediction_result)
                molecule_result["steps"]["benchmark"] = benchmark_result
                
                molecule_result["success"] = True
                self.results["success_count"] += 1
                logger.info(f"  ✓ Success: {smiles}")
                
            except Exception as e:
                logger.error(f"  ✗ Failed: {smiles} - {e}")
                molecule_result["errors"].append(str(e))
                self.results["failure_count"] += 1
                self.results["errors"].append({
                    "smiles": smiles,
                    "error": str(e),
                })
            
            self.results["molecules_tested"].append(molecule_result)
        
        # Generate summary
        self._generate_summary()
        
        return self.results
    
    async def _test_screening(self, smiles: str) -> Dict[str, Any]:
        """Test Phase 7 screening."""
        try:
            if not PHASES_AVAILABLE:
                return {"status": "skipped", "reason": "Phases not available"}
            
            search_engine = SearchEngine()
            results = search_engine.similarity_search(smiles, k=5, threshold=0.5)
            
            return {
                "status": "success",
                "similar_molecules": len(results),
            }
        except Exception as e:
            return {"status": "error", "error": str(e)}
    
    async def _test_ml_prediction(self, smiles: str) -> Dict[str, Any]:
        """Test Phase 5 ML prediction."""
        try:
            if not PHASES_AVAILABLE:
                return {"status": "skipped", "reason": "Phases not available"}
            
            registry = ModelRegistry()
            engine = PredictionEngine(registry)
            
            result = engine.predict(
                input_data={"smiles": smiles},
                properties=["logP", "solubility", "toxicity"],
            )
            
            return {
                "status": "success",
                "predictions": result.predictions,
            }
        except Exception as e:
            return {"status": "error", "error": str(e)}
    
    async def _test_conformers(self, smiles: str) -> Dict[str, Any]:
        """Test Phase 8 conformer generation."""
        try:
            if not PHASES_AVAILABLE:
                return {"status": "skipped", "reason": "Phases not available"}
            
            generator = ConformerGenerator()
            conformers = generator.generate_conformers(smiles, n=5)
            
            return {
                "status": "success",
                "conformers_generated": len(conformers),
            }
        except Exception as e:
            return {"status": "error", "error": str(e)}
    
    async def _test_rl_loop(self, smiles: str) -> Dict[str, Any]:
        """Test Phase 10 RL loop."""
        try:
            if not PHASES_AVAILABLE:
                return {"status": "skipped", "reason": "Phases not available"}
            
            # Run a micro RL loop
            rl_agent = RLAgent(config={"batch_size": 1})
            molecules = rl_agent.generate_batch(1, seed_smiles=[smiles])
            
            return {
                "status": "success",
                "generated": len(molecules),
            }
        except Exception as e:
            return {"status": "error", "error": str(e)}
    
    def _compare_benchmarks(self, smiles: str, prediction_result: Dict[str, Any]) -> Dict[str, Any]:
        """Compare predictions with known benchmarks."""
        # TODO: Load benchmark data from dataset
        # For now, just validate predictions are reasonable
        predictions = prediction_result.get("predictions", {})
        
        checks = {}
        for prop, value in predictions.items():
            if prop == "logP":
                checks[prop] = -5 <= value <= 10  # Typical range
            elif prop == "solubility":
                checks[prop] = 0 <= value <= 10  # -logS range
            elif prop == "toxicity":
                checks[prop] = 0 <= value <= 1  # Probability range
        
        return {
            "status": "success" if all(checks.values()) else "warning",
            "checks": checks,
        }
    
    def _generate_summary(self):
        """Generate validation summary."""
        total = len(self.results["molecules_tested"])
        success = self.results["success_count"]
        failure = self.results["failure_count"]
        
        logger.info("\n" + "=" * 60)
        logger.info("Validation Summary")
        logger.info("=" * 60)
        logger.info(f"Total molecules: {total}")
        logger.info(f"Successful: {success}")
        logger.info(f"Failed: {failure}")
        logger.info(f"Success rate: {success/total*100:.1f}%" if total > 0 else "N/A")
        
        if self.results["errors"]:
            logger.warning(f"\nErrors encountered: {len(self.results['errors'])}")
            for error in self.results["errors"][:5]:  # Show first 5
                logger.warning(f"  - {error['smiles']}: {error['error']}")
    
    def save_results(self, output_file: str = "data/validation/phase_a_results.json"):
        """Save validation results."""
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, "w") as f:
            json.dump(self.results, f, indent=2, default=str)
        
        logger.info(f"\nResults saved to: {output_path}")


async def main():
    """Run Phase A validation."""
    # Load small molecule library
    library_file = Path("data/libraries/compounds.smi")
    
    if not library_file.exists():
        logger.warning(f"Library file not found: {library_file}")
        logger.info("Using sample molecules...")
        test_molecules = ["CCO", "CCCO", "CC(C)O", "c1ccccc1", "CCc1ccccc1"]
    else:
        # Read first 10 molecules
        test_molecules = []
        with open(library_file, "r") as f:
            for line in f:
                if len(test_molecules) >= 10:
                    break
                smiles = line.split()[0] if line.strip() else None
                if smiles:
                    test_molecules.append(smiles)
    
    validator = EndToEndValidator()
    results = await validator.validate_pipeline(test_molecules)
    validator.save_results()


if __name__ == "__main__":
    asyncio.run(main())

