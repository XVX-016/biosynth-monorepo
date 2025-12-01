"""
Main Training Script for Phase 5 ML Models

Orchestrates the full training pipeline:
1. Load datasets
2. Configure models
3. Train models
4. Evaluate models
5. Register models
6. Test API
"""

import sys
from pathlib import Path
import logging

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

from scripts.load_datasets import load_datasets
from scripts.configure_model import configure_models
from scripts.train_models import train_models
from scripts.evaluate_models import evaluate_models
from scripts.register_models import register_models
from scripts.test_predictions import test_predictions

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main():
    """Run full training pipeline."""
    logger.info("=" * 60)
    logger.info("Phase 5 ML Model Training Pipeline")
    logger.info("=" * 60)
    
    # Step 1: Load datasets
    logger.info("\n[Step 1] Loading datasets...")
    train_data, val_data, test_data = load_datasets()
    
    # Step 2: Configure models
    logger.info("\n[Step 2] Configuring models...")
    model_configs = configure_models(train_data)
    
    # Step 3: Train models
    logger.info("\n[Step 3] Training models...")
    trained_models = train_models(model_configs, train_data, val_data)
    
    # Step 4: Evaluate models
    logger.info("\n[Step 4] Evaluating models...")
    evaluation_results = evaluate_models(trained_models, test_data)
    
    # Step 5: Register models
    logger.info("\n[Step 5] Registering models...")
    register_models(trained_models, evaluation_results)
    
    # Step 6: Test API
    logger.info("\n[Step 6] Testing API predictions...")
    test_predictions()
    
    logger.info("\n" + "=" * 60)
    logger.info("Training pipeline complete!")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()

