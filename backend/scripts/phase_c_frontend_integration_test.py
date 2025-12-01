"""
Phase C: Backend → Frontend Integration Test

Goal: Connect real outputs to UI without worrying about polish.

Tasks:
1. Feed Phase 10 RL candidates into frontend 3D viewer
2. Display Phase 5 ML predictions and attention maps
3. Show Phase 7 screening and Phase 8 conformer info panels
4. Test interactivity (click, hover, selection)
"""

import sys
from pathlib import Path
import logging
import json
import requests
from typing import Dict, Any, List

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

API_BASE = "http://localhost:8000"


class FrontendIntegrationTester:
    """Tests backend-frontend integration."""
    
    def __init__(self):
        self.results = {
            "timestamp": None,
            "api_tests": {},
            "data_flow_tests": {},
            "errors": [],
        }
    
    def test_api_endpoints(self) -> Dict[str, Any]:
        """Test all Phase 10 API endpoints."""
        logger.info("=" * 60)
        logger.info("Phase C: Backend → Frontend Integration Test")
        logger.info("=" * 60)
        
        test_smiles = "CCO"
        
        # Test 1: Generate molecules
        logger.info("\nTest 1: Generate molecules")
        try:
            response = requests.post(
                f"{API_BASE}/api/phase10/generate",
                json={"n": 5, "method": "rl"},
                timeout=10,
            )
            if response.status_code == 200:
                data = response.json()
                logger.info(f"  ✓ Generated {data.get('count', 0)} molecules")
                self.results["api_tests"]["generate"] = {"status": "success", "data": data}
            else:
                logger.error(f"  ✗ Failed: {response.status_code}")
                self.results["api_tests"]["generate"] = {"status": "error", "code": response.status_code}
        except Exception as e:
            logger.error(f"  ✗ Error: {e}")
            self.results["api_tests"]["generate"] = {"status": "error", "error": str(e)}
        
        # Test 2: Evaluate molecule
        logger.info("\nTest 2: Evaluate molecule")
        try:
            response = requests.post(
                f"{API_BASE}/api/phase10/evaluate",
                json={
                    "smiles": test_smiles,
                    "compute_ml": True,
                    "compute_screening": True,
                },
                timeout=10,
            )
            if response.status_code == 200:
                data = response.json()
                logger.info(f"  ✓ Evaluation successful, reward: {data.get('reward', 0):.4f}")
                self.results["api_tests"]["evaluate"] = {"status": "success", "data": data}
            else:
                logger.error(f"  ✗ Failed: {response.status_code}")
                self.results["api_tests"]["evaluate"] = {"status": "error", "code": response.status_code}
        except Exception as e:
            logger.error(f"  ✗ Error: {e}")
            self.results["api_tests"]["evaluate"] = {"status": "error", "error": str(e)}
        
        # Test 3: Get top candidates
        logger.info("\nTest 3: Get top candidates")
        try:
            response = requests.get(f"{API_BASE}/api/phase10/top_candidates?n=5", timeout=10)
            if response.status_code == 200:
                data = response.json()
                logger.info(f"  ✓ Retrieved {len(data.get('candidates', []))} candidates")
                self.results["api_tests"]["top_candidates"] = {"status": "success", "data": data}
            else:
                logger.error(f"  ✗ Failed: {response.status_code}")
                self.results["api_tests"]["top_candidates"] = {"status": "error", "code": response.status_code}
        except Exception as e:
            logger.error(f"  ✗ Error: {e}")
            self.results["api_tests"]["top_candidates"] = {"status": "error", "error": str(e)}
        
        # Test 4: Get statistics
        logger.info("\nTest 4: Get statistics")
        try:
            response = requests.get(f"{API_BASE}/api/phase10/statistics", timeout=10)
            if response.status_code == 200:
                data = response.json()
                logger.info(f"  ✓ Statistics retrieved")
                self.results["api_tests"]["statistics"] = {"status": "success", "data": data}
            else:
                logger.error(f"  ✗ Failed: {response.status_code}")
                self.results["api_tests"]["statistics"] = {"status": "error", "code": response.status_code}
        except Exception as e:
            logger.error(f"  ✗ Error: {e}")
            self.results["api_tests"]["statistics"] = {"status": "error", "error": str(e)}
        
        return self.results
    
    def test_data_flow(self) -> Dict[str, Any]:
        """Test data flow through pipeline."""
        logger.info("\n" + "=" * 60)
        logger.info("Data Flow Test")
        logger.info("=" * 60)
        
        test_smiles = "CCO"
        
        # Test: SMILES → ML Prediction → Conformers → RL Loop
        logger.info(f"\nTesting data flow for: {test_smiles}")
        
        flow_steps = {}
        
        # Step 1: ML Prediction
        logger.info("  → ML Prediction")
        try:
            response = requests.post(
                f"{API_BASE}/api/predict/property",
                json={"smiles": test_smiles, "properties": ["logP", "solubility"]},
                timeout=10,
            )
            if response.status_code == 200:
                flow_steps["ml_prediction"] = {"status": "success", "data": response.json()}
                logger.info("    ✓ ML prediction successful")
            else:
                flow_steps["ml_prediction"] = {"status": "error", "code": response.status_code}
        except Exception as e:
            flow_steps["ml_prediction"] = {"status": "error", "error": str(e)}
        
        # Step 2: Conformers
        logger.info("  → Conformers")
        try:
            response = requests.post(
                f"{API_BASE}/api/conformers/generate",
                json={"smiles": test_smiles, "n": 3},
                timeout=10,
            )
            if response.status_code == 200:
                flow_steps["conformers"] = {"status": "success", "data": response.json()}
                logger.info("    ✓ Conformers generated")
            else:
                flow_steps["conformers"] = {"status": "error", "code": response.status_code}
        except Exception as e:
            flow_steps["conformers"] = {"status": "error", "error": str(e)}
        
        # Step 3: RL Loop
        logger.info("  → RL Loop")
        try:
            response = requests.post(
                f"{API_BASE}/api/phase10/run_loop",
                json={"max_iterations": 2, "batch_size": 3},
                timeout=30,
            )
            if response.status_code == 200:
                flow_steps["rl_loop"] = {"status": "success", "data": response.json()}
                logger.info("    ✓ RL loop executed")
            else:
                flow_steps["rl_loop"] = {"status": "error", "code": response.status_code}
        except Exception as e:
            flow_steps["rl_loop"] = {"status": "error", "error": str(e)}
        
        self.results["data_flow_tests"] = flow_steps
        
        # Summary
        success_count = sum(1 for step in flow_steps.values() if step.get("status") == "success")
        logger.info(f"\nData flow summary: {success_count}/{len(flow_steps)} steps successful")
        
        return flow_steps
    
    def save_results(self, output_file: str = "data/validation/phase_c_results.json"):
        """Save integration test results."""
        from datetime import datetime
        self.results["timestamp"] = datetime.now().isoformat()
        
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, "w") as f:
            json.dump(self.results, f, indent=2, default=str)
        
        logger.info(f"\nResults saved to: {output_path}")


def main():
    """Run Phase C integration tests."""
    tester = FrontendIntegrationTester()
    
    # Test API endpoints
    tester.test_api_endpoints()
    
    # Test data flow
    tester.test_data_flow()
    
    # Save results
    tester.save_results()
    
    logger.info("\n" + "=" * 60)
    logger.info("Integration test complete!")
    logger.info("=" * 60)
    logger.info("\nNote: Frontend UI testing requires manual verification:")
    logger.info("  1. Start backend: python -m uvicorn app:app")
    logger.info("  2. Start frontend: cd frontend && npm run dev")
    logger.info("  3. Navigate to /phase10")
    logger.info("  4. Test interactivity manually")


if __name__ == "__main__":
    main()

