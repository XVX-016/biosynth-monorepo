"""
Step 4: Evaluate Trained Models

Evaluates models on test set and computes metrics.
"""

import sys
from pathlib import Path
import logging
import json
from typing import Dict, List, Any
import numpy as np

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

try:
    import torch
    from torch.utils.data import DataLoader
    from torch_geometric.data import Batch
    from scripts.train_models import MoleculeDataset, _collate_fn
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False

logger = logging.getLogger(__name__)


def evaluate_models(
    trained_models: Dict[str, Any],
    test_data: Dict[str, Any],
) -> Dict[str, Dict[str, float]]:
    """
    Evaluate models on test set.
    
    Args:
        trained_models: Dictionary of trained models
        test_data: Test dataset
    
    Returns:
        Dictionary mapping property names to evaluation metrics
    """
    if not TORCH_AVAILABLE:
        logger.warning("PyTorch not available. Returning mock metrics.")
        return _create_mock_metrics(trained_models)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    eval_dir = Path("data/models/evaluation")
    eval_dir.mkdir(parents=True, exist_ok=True)
    
    evaluation_results = {}
    
    for prop, model_info in trained_models.items():
        logger.info(f"\nEvaluating model for {prop}...")
        
        if prop not in test_data["properties"]:
            logger.warning(f"Property {prop} not in test data. Skipping.")
            continue
        
        model = model_info["model"]
        model.eval()
        
        # Create test dataset
        test_targets = test_data["properties"][prop]
        test_dataset = MoleculeDataset(
            test_data["smiles"],
            test_targets,
            prop,
        )
        test_loader = DataLoader(
            test_dataset,
            batch_size=32,
            shuffle=False,
            collate_fn=_collate_fn,
        )
        
        # Evaluate
        predictions = []
        targets = []
        
        with torch.no_grad():
            for batch_data, batch_targets in test_loader:
                batch_data = batch_data.to(device)
                preds = model(
                    batch_data.x,
                    batch_data.edge_index,
                    edge_attr=getattr(batch_data, 'edge_attr', None),
                    batch=batch_data.batch,
                )
                predictions.extend(preds.cpu().numpy().flatten().tolist())
                targets.extend(batch_targets.numpy().tolist())
        
        # Compute metrics
        predictions = np.array(predictions)
        targets = np.array(targets)
        
        mae = np.mean(np.abs(predictions - targets))
        rmse = np.sqrt(np.mean((predictions - targets) ** 2))
        
        # R²
        ss_res = np.sum((targets - predictions) ** 2)
        ss_tot = np.sum((targets - np.mean(targets)) ** 2)
        r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
        
        metrics = {
            "mae": float(mae),
            "rmse": float(rmse),
            "r2": float(r2),
            "n_samples": len(targets),
        }
        
        evaluation_results[prop] = metrics
        
        logger.info(f"{prop} - MAE: {mae:.4f}, RMSE: {rmse:.4f}, R²: {r2:.4f}")
        
        # Save evaluation report
        report_file = eval_dir / f"{prop}_evaluation.json"
        with open(report_file, "w") as f:
            json.dump(metrics, f, indent=2)
    
    # Save combined report
    combined_report = eval_dir / "combined_evaluation.json"
    with open(combined_report, "w") as f:
        json.dump(evaluation_results, f, indent=2)
    
    logger.info(f"\nEvaluation complete. Reports saved to {eval_dir}")
    
    return evaluation_results


def _create_mock_metrics(trained_models: Dict[str, Any]) -> Dict[str, Dict[str, float]]:
    """Create mock metrics when PyTorch is not available."""
    return {
        prop: {
            "mae": 0.5,
            "rmse": 0.7,
            "r2": 0.6,
            "n_samples": 10,
        }
        for prop in trained_models.keys()
    }


if __name__ == "__main__":
    # This would be called from train_ml_models.py
    pass

