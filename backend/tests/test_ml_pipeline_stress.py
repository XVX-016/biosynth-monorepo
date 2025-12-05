"""
ML Pipeline Stress Tests

Phase 16: Full Regression Test & Final Cleanup

Tests ML prediction pipeline with:
- Large batch predictions (1000 molecules)
- Various molecule types
- Performance profiling
"""

import pytest
import requests
import time
import logging
from typing import List, Dict, Any

logger = logging.getLogger(__name__)

API_BASE = "http://localhost:8000"

# Generate test SMILES for stress testing
def generate_test_smiles(count: int = 100) -> List[str]:
    """Generate a list of test SMILES."""
    base_smiles = [
        "C", "CC", "CCC", "CCCC", "CCCCC",  # Alkanes
        "C=C", "CC=CC", "C=CCC",  # Alkenes
        "C#C", "CC#CC",  # Alkynes
        "CCO", "CCCO", "CCOC",  # Alcohols and ethers
        "c1ccccc1", "c1ccccc1C",  # Aromatics
        "CC(=O)O", "CC(=O)CC",  # Ketones
        "CCOC(=O)O",  # Esters
    ]
    
    # Repeat and vary
    smiles_list = []
    for i in range(count):
        smiles_list.append(base_smiles[i % len(base_smiles)])
    
    return smiles_list


class TestMLPipelineStress:
    """Stress tests for ML prediction pipeline."""

    @pytest.fixture(autouse=True)
    def setup(self):
        """Setup for each test."""
        self.base_url = API_BASE

    def test_batch_prediction_100(self):
        """Test batch prediction with 100 molecules."""
        smiles_list = generate_test_smiles(100)
        
        inputs = [{"smiles": smiles} for smiles in smiles_list]
        
        start_time = time.time()
        response = requests.post(
            f"{self.base_url}/api/predict/batch",
            json={"inputs": inputs},
            timeout=120,  # 2 minutes for 100 molecules
        )
        elapsed = time.time() - start_time
        
        assert response.status_code == 200, f"Batch prediction failed: {response.status_code}"
        data = response.json()
        assert "results" in data, "Missing results in response"
        assert len(data["results"]) == 100, f"Expected 100 results, got {len(data['results'])}"
        
        logger.info(f"✓ Batch prediction (100 molecules): {elapsed:.2f}s ({elapsed/100*1000:.2f}ms per molecule)")

    def test_batch_prediction_1000(self):
        """Test batch prediction with 1000 molecules."""
        smiles_list = generate_test_smiles(1000)
        
        inputs = [{"smiles": smiles} for smiles in smiles_list]
        
        start_time = time.time()
        response = requests.post(
            f"{self.base_url}/api/predict/batch",
            json={"inputs": inputs},
            timeout=600,  # 10 minutes for 1000 molecules
        )
        elapsed = time.time() - start_time
        
        assert response.status_code == 200, f"Batch prediction failed: {response.status_code}"
        data = response.json()
        assert "results" in data, "Missing results in response"
        assert len(data["results"]) == 1000, f"Expected 1000 results, got {len(data['results'])}"
        
        logger.info(f"✓ Batch prediction (1000 molecules): {elapsed:.2f}s ({elapsed/1000*1000:.2f}ms per molecule)")

    def test_property_prediction_performance(self):
        """Test property prediction performance."""
        test_smiles = ["CCO", "c1ccccc1", "CC(=O)O", "CCOC", "CCCC"]
        
        times = []
        for smiles in test_smiles:
            start_time = time.time()
            response = requests.post(
                f"{self.base_url}/api/predict/property",
                json={"smiles": smiles, "properties": ["logP", "solubility", "toxicity"]},
                timeout=30,
            )
            elapsed = time.time() - start_time
            times.append(elapsed)
            
            assert response.status_code == 200, f"Prediction failed for {smiles}"
            data = response.json()
            assert "predictions" in data, "Missing predictions in response"
        
        avg_time = sum(times) / len(times)
        max_time = max(times)
        min_time = min(times)
        
        logger.info(f"✓ Property prediction performance:")
        logger.info(f"  Average: {avg_time*1000:.2f}ms")
        logger.info(f"  Min: {min_time*1000:.2f}ms")
        logger.info(f"  Max: {max_time*1000:.2f}ms")
        
        # Performance assertion (should be < 5s per prediction)
        assert avg_time < 5.0, f"Average prediction time too high: {avg_time:.2f}s"

    def test_attention_map_performance(self):
        """Test attention map generation performance."""
        test_smiles = ["CCO", "c1ccccc1", "CC(=O)O"]
        
        times = []
        for smiles in test_smiles:
            start_time = time.time()
            response = requests.post(
                f"{self.base_url}/api/predict/attention-map",
                json={"smiles": smiles},
                timeout=30,
            )
            elapsed = time.time() - start_time
            times.append(elapsed)
            
            assert response.status_code == 200, f"Attention map failed for {smiles}"
            data = response.json()
            assert "attention_map" in data or "attentions" in data, "Missing attention data"
        
        avg_time = sum(times) / len(times)
        logger.info(f"✓ Attention map performance: {avg_time*1000:.2f}ms average")
        
        # Performance assertion
        assert avg_time < 10.0, f"Average attention map time too high: {avg_time:.2f}s"

    def test_concurrent_requests(self):
        """Test handling of concurrent requests."""
        import concurrent.futures
        
        smiles_list = ["CCO", "c1ccccc1", "CC(=O)O", "CCOC", "CCCC"]
        
        def make_request(smiles: str):
            response = requests.post(
                f"{self.base_url}/api/predict/property",
                json={"smiles": smiles},
                timeout=30,
            )
            return response.status_code == 200
        
        start_time = time.time()
        with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
            results = list(executor.map(make_request, smiles_list))
        elapsed = time.time() - start_time
        
        assert all(results), "Some concurrent requests failed"
        logger.info(f"✓ Concurrent requests (5): {elapsed:.2f}s")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

