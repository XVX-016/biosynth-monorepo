"""
Step 6: Test API Predictions

Tests prediction endpoints with real models.
"""

import sys
from pathlib import Path
import logging
import requests
import json

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

logger = logging.getLogger(__name__)

API_BASE = "http://localhost:8000"


def test_predictions():
    """
    Test prediction API endpoints.
    
    Tests:
    - POST /api/predict/property
    - Confirms output values are realistic
    - Confirms attention maps are generated (if applicable)
    """
    logger.info("Testing prediction API endpoints...")
    
    test_smiles = "CCO"  # Ethanol
    test_properties = ["logP", "toxicity"]
    
    # Test 1: Property prediction
    logger.info(f"\nTest 1: Property prediction for {test_smiles}")
    try:
        response = requests.post(
            f"{API_BASE}/api/predict/property",
            json={
                "smiles": test_smiles,
                "properties": test_properties,
            },
            timeout=10,
        )
        
        if response.status_code == 200:
            result = response.json()
            logger.info(f"✓ Prediction successful")
            logger.info(f"  Predictions: {result.get('predictions', {})}")
            
            # Check if values are realistic
            predictions = result.get("predictions", {})
            for prop, value in predictions.items():
                if prop == "logP":
                    if -5 <= value <= 10:
                        logger.info(f"  ✓ {prop} value is realistic: {value:.4f}")
                    else:
                        logger.warning(f"  ⚠ {prop} value may be unrealistic: {value:.4f}")
                elif prop == "toxicity":
                    if 0 <= value <= 1:
                        logger.info(f"  ✓ {prop} value is realistic: {value:.4f}")
                    else:
                        logger.warning(f"  ⚠ {prop} value may be unrealistic: {value:.4f}")
        else:
            logger.error(f"✗ Prediction failed: {response.status_code}")
            logger.error(f"  Response: {response.text}")
    except requests.exceptions.ConnectionError:
        logger.warning("✗ API server not running. Start backend with: python -m uvicorn app:app")
    except Exception as e:
        logger.error(f"✗ Test failed: {e}")
    
    # Test 2: Attention maps
    logger.info(f"\nTest 2: Attention map generation")
    try:
        response = requests.post(
            f"{API_BASE}/api/predict/attention-map",
            json={
                "smiles": test_smiles,
                "properties": test_properties,
            },
            timeout=10,
        )
        
        if response.status_code == 200:
            result = response.json()
            logger.info(f"✓ Attention map generated")
            
            if "attentions" in result:
                logger.info(f"  ✓ Attention weights present")
            else:
                logger.warning(f"  ⚠ No attention weights in response")
        else:
            logger.warning(f"✗ Attention map endpoint failed: {response.status_code}")
    except requests.exceptions.ConnectionError:
        logger.warning("✗ API server not running")
    except Exception as e:
        logger.warning(f"✗ Attention test failed: {e}")
    
    logger.info("\nAPI testing complete!")


if __name__ == "__main__":
    test_predictions()

