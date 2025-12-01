"""
Direct Model Prediction Testing

Tests trained models directly without requiring API server.
"""

import sys
from pathlib import Path
import logging

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def test_predictions():
    """Test model predictions directly."""
    logger.info("=" * 60)
    logger.info("Direct Model Prediction Testing")
    logger.info("=" * 60)
    
    try:
        from ml.prediction_engine import PredictionEngine
        from ml.registry import ModelRegistry
        
        logger.info("\n✓ Prediction engine imported")
        
        # Initialize with correct registry path
        registry = ModelRegistry(registry_path="data/models/registry.json")
        engine = PredictionEngine(registry)
        
        # List available models
        models = registry.list_models()
        logger.info(f"\nAvailable models: {len(models)}")
        for model in models:
            logger.info(f"  - {model['id']}: {model.get('description', 'N/A')}")
        
        # Test molecules
        test_molecules = {
            "CCO": "Ethanol",
            "CCCO": "Propanol",
            "CC(C)O": "Isopropanol",
            "c1ccccc1": "Benzene",
            "CCc1ccccc1": "Ethylbenzene",
        }
        
        logger.info(f"\nTesting predictions for {len(test_molecules)} molecules...")
        
        results = []
        for smiles, name in test_molecules.items():
            try:
                logger.info(f"\n  Testing: {name} ({smiles})")
                
                # Test with all properties
                result = engine.predict(
                    input_data={"smiles": smiles},
                    properties=["logP", "solubility", "toxicity"],
                    return_attention=False,  # Skip attention for faster testing
                )
                
                # Also test individual property predictions
                logp_result = engine.predict(
                    input_data={"smiles": smiles},
                    properties=["logP"],
                    model_id="attention-gnn-logp",
                )
                
                logger.info(f"    Predictions:")
                for prop, value in result.predictions.items():
                    logger.info(f"      {prop}: {value:.4f}")
                
                results.append({
                    "smiles": smiles,
                    "name": name,
                    "predictions": result.predictions,
                    "success": True,
                })
                
            except Exception as e:
                logger.error(f"    ✗ Failed: {e}")
                results.append({
                    "smiles": smiles,
                    "name": name,
                    "error": str(e),
                    "success": False,
                })
        
        # Summary
        success_count = sum(1 for r in results if r.get("success", False))
        logger.info("\n" + "=" * 60)
        logger.info("Prediction Test Summary")
        logger.info("=" * 60)
        logger.info(f"Total molecules: {len(results)}")
        logger.info(f"Successful: {success_count}")
        logger.info(f"Failed: {len(results) - success_count}")
        
        if success_count > 0:
            logger.info("\n✓ Model predictions working!")
        
        return results
        
    except ImportError as e:
        logger.error(f"✗ Import failed: {e}")
        logger.info("\nNote: PyTorch/PyG may not be available")
        return []
    except Exception as e:
        logger.error(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return []


if __name__ == "__main__":
    test_predictions()

