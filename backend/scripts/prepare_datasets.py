"""
Dataset Preparation Script

Main entry point for preparing datasets for ML training and Phase 10 RL.
Runs all preparation steps in sequence.
"""

import sys
from pathlib import Path

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

from scripts.collect_property_datasets import collect_property_datasets
from scripts.prepare_screening_libraries import prepare_screening_libraries
from scripts.standardize_molecules import standardize_molecules
from scripts.split_datasets import split_datasets
from scripts.generate_fingerprints import generate_fingerprints
from scripts.generate_mock_rewards import generate_mock_rewards

import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main():
    """Run all dataset preparation steps."""
    logger.info("Starting dataset preparation...")
    
    # Step 1: Collect property datasets
    logger.info("Step 1: Collecting property datasets...")
    collect_property_datasets()
    
    # Step 2: Prepare screening libraries
    logger.info("Step 2: Preparing screening libraries...")
    prepare_screening_libraries()
    
    # Step 3: Standardize molecules
    logger.info("Step 3: Standardizing molecules...")
    standardize_molecules()
    
    # Step 4: Split datasets
    logger.info("Step 4: Splitting datasets...")
    split_datasets()
    
    # Step 5: Generate fingerprints
    logger.info("Step 5: Generating fingerprints...")
    generate_fingerprints()
    
    # Step 6: Generate mock rewards (optional)
    logger.info("Step 6: Generating mock rewards...")
    generate_mock_rewards()
    
    logger.info("Dataset preparation complete!")


if __name__ == "__main__":
    main()

