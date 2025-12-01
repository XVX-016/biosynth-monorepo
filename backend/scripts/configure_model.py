"""
Step 2: Configure Attention-GNN Models

Initializes attention-GNN models for each property.
"""

import sys
from pathlib import Path
import logging
from typing import Dict, List, Any

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

try:
    from ml.gat_model import create_model
    from ml.featurize import featurize_smiles
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False
    logging.warning("PyTorch/PyG not available. Models will be mock.")

logger = logging.getLogger(__name__)


def configure_models(train_data: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    """
    Configure models for each property.
    
    Args:
        train_data: Training dataset
    
    Returns:
        Dictionary mapping property names to model configurations
    """
    logger.info("Configuring Attention-GNN models...")
    
    # Determine input dimensions from sample molecule
    node_feat_dim = 64  # Default
    edge_feat_dim = 8  # Default
    
    if TORCH_AVAILABLE and train_data.get("smiles"):
        try:
            # Featurize a sample molecule to get dimensions
            sample_smiles = train_data["smiles"][0]
            result = featurize_smiles(sample_smiles)
            # featurize_smiles returns (data, node_mapping) tuple
            if isinstance(result, tuple):
                data, _ = result
            else:
                data = result
            if data and hasattr(data, 'x') and data.x is not None:
                node_feat_dim = data.x.shape[1]
                logger.info(f"Detected node_feat_dim: {node_feat_dim}")
            if data and hasattr(data, 'edge_attr') and data.edge_attr is not None and len(data.edge_attr) > 0:
                edge_feat_dim = data.edge_attr.shape[1]
                logger.info(f"Detected edge_feat_dim: {edge_feat_dim}")
        except Exception as e:
            logger.warning(f"Could not determine dimensions from sample: {e}")
            # Use actual dimensions from featurizer
            node_feat_dim = 12  # From compute_node_features
            edge_feat_dim = 8   # From compute_edge_features
    
    # Properties to train
    properties = ["logP", "solubility", "toxicity"]
    
    model_configs = {}
    for prop in properties:
        model_configs[prop] = {
            "property": prop,
            "node_feat_dim": node_feat_dim,
            "edge_feat_dim": edge_feat_dim,
            "hidden_dim": 128,
            "out_dim": 1,  # Single property per model
            "num_layers": 3,
            "heads": 4,
            "dropout": 0.1,
            "learning_rate": 0.001,
            "batch_size": 32,
            "epochs": 50,
        }
        logger.info(f"Configured model for {prop}: {node_feat_dim}D nodes, {edge_feat_dim}D edges")
    
    return model_configs


if __name__ == "__main__":
    # Mock data for testing
    mock_data = {
        "smiles": ["CCO", "CCCO"],
        "properties": {"logP": [0.0, 0.5]},
    }
    configs = configure_models(mock_data)
    print(configs)

