"""
Step 6: Generate Mock Rewards for RL Testing

Creates temporary CSV with SMILES + random reward values
for Phase 10 RL loop testing before real ML models are trained.
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

try:
    from backend.chem.utils.validators import validate_smiles
except ImportError:
    # Fallback: simple validation
    def validate_smiles(smiles: str) -> bool:
        return bool(smiles and len(smiles) > 0)

logger = logging.getLogger(__name__)


def generate_mock_rewards():
    """
    Generate mock reward scores for molecules.
    
    Creates a CSV with SMILES and random reward values for
    testing Phase 10 RL loop before real ML models are trained.
    """
    data_dir = Path("data/datasets")
    cleaned_dir = data_dir / "cleaned"
    
    # Try to get molecules from cleaned dataset
    input_file = cleaned_dir / "properties_cleaned.csv"
    if not input_file.exists():
        input_file = data_dir / "properties.csv"
    
    output_file = data_dir / "mock_rewards.csv"
    
    if not input_file.exists():
        logger.warning(f"Input file not found: {input_file}")
        logger.info("Creating mock rewards from sample molecules...")
        _create_sample_mock_rewards(output_file)
        return
    
    logger.info(f"Generating mock rewards from {input_file}")
    
    random.seed(42)  # Reproducible
    
    mock_rewards = []
    with open(input_file, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            smiles = row.get("smiles", "").strip()
            if not smiles or not validate_smiles(smiles):
                continue
            
            # Generate mock reward (0.0 to 1.0)
            # Slightly bias towards higher rewards for some molecules
            base_reward = random.uniform(0.0, 1.0)
            
            # Add some structure-based bias (longer molecules get slightly lower rewards)
            if len(smiles) > 20:
                base_reward *= 0.8
            
            mock_rewards.append({
                "smiles": smiles,
                "reward": round(base_reward, 4),
            })
    
    # Write mock rewards
    if mock_rewards:
        fieldnames = ["smiles", "reward"]
        with open(output_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(mock_rewards)
        
        logger.info(f"Generated {len(mock_rewards)} mock rewards: {output_file}")
    else:
        logger.error("No mock rewards generated!")


def _create_sample_mock_rewards(output_file: Path):
    """Create sample mock rewards if no input file exists."""
    sample_molecules = [
        "CCO",
        "CC(=O)O",
        "c1ccccc1",
        "CCc1ccccc1",
        "CC(C)O",
        "CCN(CC)CC",
        "CC(=O)OC",
        "CCCCCCCCCC",
    ]
    
    random.seed(42)
    mock_rewards = []
    for smiles in sample_molecules:
        mock_rewards.append({
            "smiles": smiles,
            "reward": round(random.uniform(0.0, 1.0), 4),
        })
    
    fieldnames = ["smiles", "reward"]
    with open(output_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(mock_rewards)
    
    logger.info(f"Created sample mock rewards: {output_file} ({len(mock_rewards)} molecules)")


if __name__ == "__main__":
    generate_mock_rewards()

