"""
Step 4: Split Datasets

Splits dataset into training (70%), validation (15%), and test (15%) sets.
"""

import sys
from pathlib import Path
import csv
import random
import logging

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

logger = logging.getLogger(__name__)


def split_datasets():
    """
    Split dataset into train/validation/test sets.
    
    Splits:
    - Training: 70%
    - Validation: 15%
    - Test: 15%
    """
    data_dir = Path("data/datasets")
    cleaned_dir = data_dir / "cleaned"
    input_file = cleaned_dir / "properties_cleaned.csv"
    
    # Fallback to original if cleaned doesn't exist
    if not input_file.exists():
        input_file = data_dir / "properties.csv"
    
    if not input_file.exists():
        logger.warning(f"Input file not found: {input_file}")
        return
    
    logger.info(f"Splitting dataset from {input_file}")
    
    # Read all data
    rows = []
    with open(input_file, "r") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    
    if not rows:
        logger.error("No data to split!")
        return
    
    # Shuffle
    random.seed(42)  # Reproducible splits
    random.shuffle(rows)
    
    # Calculate split indices
    total = len(rows)
    train_end = int(total * 0.70)
    val_end = train_end + int(total * 0.15)
    
    train_rows = rows[:train_end]
    val_rows = rows[train_end:val_end]
    test_rows = rows[val_end:]
    
    # Write splits
    fieldnames = list(rows[0].keys())
    
    train_file = cleaned_dir / "train.csv"
    with open(train_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(train_rows)
    logger.info(f"Training set: {train_file} ({len(train_rows)} molecules)")
    
    val_file = cleaned_dir / "validation.csv"
    with open(val_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(val_rows)
    logger.info(f"Validation set: {val_file} ({len(val_rows)} molecules)")
    
    test_file = cleaned_dir / "test.csv"
    with open(test_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(test_rows)
    logger.info(f"Test set: {test_file} ({len(test_rows)} molecules)")


if __name__ == "__main__":
    split_datasets()

