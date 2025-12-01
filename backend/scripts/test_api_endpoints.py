"""
Test API Endpoints

Tests all ML prediction API endpoints with the running server.
"""

import requests
import json
import time
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

API_BASE = "http://localhost:8000"


def wait_for_server(max_retries=10, delay=2):
    """Wait for server to be ready."""
    for i in range(max_retries):
        try:
            response = requests.get(f"{API_BASE}/health", timeout=2)
            if response.status_code == 200:
                logger.info("✓ Server is ready")
                return True
        except requests.exceptions.ConnectionError:
            logger.info(f"Waiting for server... ({i+1}/{max_retries})")
            time.sleep(delay)
    return False


def test_property_prediction():
    """Test property prediction endpoint."""
    logger.info("\n" + "=" * 60)
    logger.info("Test 1: Property Prediction")
    logger.info("=" * 60)
    
    test_cases = [
        {"smiles": "CCO", "name": "Ethanol"},
        {"smiles": "c1ccccc1", "name": "Benzene"},
        {"smiles": "CCCO", "name": "Propanol"},
    ]
    
    for test in test_cases:
        try:
            response = requests.post(
                f"{API_BASE}/api/predict/property",
                json={
                    "smiles": test["smiles"],
                    "properties": ["logP", "solubility", "toxicity"],
                },
                timeout=10,
            )
            
            if response.status_code == 200:
                result = response.json()
                logger.info(f"\n✓ {test['name']} ({test['smiles']})")
                logger.info(f"  Predictions: {result.get('predictions', {})}")
            else:
                logger.error(f"✗ {test['name']}: {response.status_code} - {response.text}")
        except Exception as e:
            logger.error(f"✗ {test['name']}: {e}")


def test_attention_map():
    """Test attention map endpoint."""
    logger.info("\n" + "=" * 60)
    logger.info("Test 2: Attention Map")
    logger.info("=" * 60)
    
    try:
        response = requests.post(
            f"{API_BASE}/api/predict/attention-map",
            json={
                "smiles": "CCO",
                "properties": ["logP"],
            },
            timeout=10,
        )
        
        if response.status_code == 200:
            result = response.json()
            logger.info("✓ Attention map generated")
            if "attentions" in result:
                logger.info(f"  Attention layers: {len(result['attentions'])}")
            if "node_importance" in result:
                logger.info(f"  Node importance: {len(result['node_importance'])} nodes")
        else:
            logger.warning(f"✗ Attention map: {response.status_code} - {response.text}")
    except Exception as e:
        logger.warning(f"✗ Attention map: {e}")


def test_batch_prediction():
    """Test batch prediction endpoint."""
    logger.info("\n" + "=" * 60)
    logger.info("Test 3: Batch Prediction")
    logger.info("=" * 60)
    
    try:
        response = requests.post(
            f"{API_BASE}/api/predict/batch",
            json={
                "inputs": [
                    {"smiles": "CCO"},
                    {"smiles": "CCCO"},
                    {"smiles": "c1ccccc1"},
                ],
                "properties": ["logP"],
            },
            timeout=15,
        )
        
        if response.status_code == 200:
            result = response.json()
            logger.info(f"✓ Batch prediction: {len(result.get('results', []))} molecules")
        else:
            logger.warning(f"✗ Batch prediction: {response.status_code} - {response.text}")
    except Exception as e:
        logger.warning(f"✗ Batch prediction: {e}")


def main():
    """Run all API tests."""
    logger.info("=" * 60)
    logger.info("API Endpoint Testing")
    logger.info("=" * 60)
    
    if not wait_for_server():
        logger.error("✗ Server not available. Start with: python -m uvicorn app:app")
        return
    
    # Run tests
    test_property_prediction()
    test_attention_map()
    test_batch_prediction()
    
    logger.info("\n" + "=" * 60)
    logger.info("API Testing Complete!")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()

