"""
Step 3: Train Attention-GNN Models

Trains models for each property using training and validation sets.
"""

import sys
from pathlib import Path
import logging
from typing import Dict, List, Any, Optional
from pathlib import Path

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

try:
    import torch
    import torch.nn as nn
    from torch.utils.data import Dataset, DataLoader
    from torch_geometric.data import Data, Batch
    from ml.gat_model import create_model
    from ml.featurize import featurize_smiles
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False
    logging.warning("PyTorch/PyG not available. Training will be mock.")

logger = logging.getLogger(__name__)


class MoleculeDataset(Dataset):
    """Dataset for molecular property prediction."""
    
    def __init__(self, smiles: List[str], targets: List[float], property_name: str):
        self.smiles = smiles
        self.targets = targets
        self.property_name = property_name
    
    def __len__(self):
        return len(self.smiles)
    
    def __getitem__(self, idx):
        smiles = self.smiles[idx]
        target = self.targets[idx]
        
        # Featurize molecule
        if TORCH_AVAILABLE:
            try:
                data = featurize_smiles(smiles)
                if data is None:
                    # Return dummy data
                    data = Data(
                        x=torch.zeros(1, 64),
                        edge_index=torch.zeros(2, 0, dtype=torch.long),
                        edge_attr=torch.zeros(0, 8),
                    )
            except Exception as e:
                logger.warning(f"Featurization failed for {smiles}: {e}")
                data = Data(
                    x=torch.zeros(1, 64),
                    edge_index=torch.zeros(2, 0, dtype=torch.long),
                    edge_attr=torch.zeros(0, 8),
                )
        else:
            data = None
        
        return data, torch.tensor(target, dtype=torch.float)


def train_models(
    model_configs: Dict[str, Dict[str, Any]],
    train_data: Dict[str, Any],
    val_data: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Train models for each property.
    
    Args:
        model_configs: Model configurations
        train_data: Training dataset
        val_data: Validation dataset
    
    Returns:
        Dictionary mapping property names to trained models and metadata
    """
    if not TORCH_AVAILABLE:
        logger.warning("PyTorch not available. Returning mock models.")
        return _create_mock_models(model_configs)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logger.info(f"Using device: {device}")
    
    models_dir = Path("data/models")
    models_dir.mkdir(parents=True, exist_ok=True)
    
    trained_models = {}
    
    for prop, config in model_configs.items():
        logger.info(f"\nTraining model for {prop}...")
        
        # Get targets
        if prop not in train_data["properties"]:
            logger.warning(f"Property {prop} not in training data. Skipping.")
            continue
        
        train_targets = train_data["properties"][prop]
        val_targets = val_data["properties"].get(prop, [])
        
        # Create datasets
        train_dataset = MoleculeDataset(
            train_data["smiles"],
            train_targets,
            prop,
        )
        val_dataset = MoleculeDataset(
            val_data["smiles"],
            val_targets,
            prop,
        ) if val_targets else None
        
        # Create model
        model = create_model(
            node_feat_dim=config["node_feat_dim"],
            edge_feat_dim=config["edge_feat_dim"],
            hidden_dim=config["hidden_dim"],
            out_dim=config["out_dim"],
            num_layers=config["num_layers"],
            heads=config["heads"],
        ).to(device)
        
        # Training setup
        optimizer = torch.optim.Adam(model.parameters(), lr=config["learning_rate"])
        criterion = nn.MSELoss()
        
        # Data loaders
        train_loader = DataLoader(
            train_dataset,
            batch_size=config["batch_size"],
            shuffle=True,
            collate_fn=_collate_fn,
        )
        val_loader = DataLoader(
            val_dataset,
            batch_size=config["batch_size"],
            shuffle=False,
            collate_fn=_collate_fn,
        ) if val_dataset else None
        
        # Training loop
        best_val_loss = float('inf')
        best_model_path = models_dir / f"{prop}_best.pt"
        
        for epoch in range(config["epochs"]):
            # Train
            model.train()
            train_loss = 0.0
            for batch_data, batch_targets in train_loader:
                batch_data = batch_data.to(device)
                batch_targets = batch_targets.to(device)
                
                optimizer.zero_grad()
                preds = model(
                    batch_data.x,
                    batch_data.edge_index,
                    edge_attr=getattr(batch_data, 'edge_attr', None),
                    batch=batch_data.batch,
                )
                loss = criterion(preds.view(-1), batch_targets)
                loss.backward()
                optimizer.step()
                train_loss += loss.item()
            
            train_loss /= len(train_loader)
            
            # Validate
            if val_loader:
                model.eval()
                val_loss = 0.0
                with torch.no_grad():
                    for batch_data, batch_targets in val_loader:
                        batch_data = batch_data.to(device)
                        batch_targets = batch_targets.to(device)
                        preds = model(
                            batch_data.x,
                            batch_data.edge_index,
                            edge_attr=getattr(batch_data, 'edge_attr', None),
                            batch=batch_data.batch,
                        )
                        loss = criterion(preds.view(-1), batch_targets)
                        val_loss += loss.item()
                
                val_loss /= len(val_loader)
                
                # Save best model
                if val_loss < best_val_loss:
                    best_val_loss = val_loss
                    torch.save(model.state_dict(), best_model_path)
                    logger.info(f"Epoch {epoch+1}: Saved best model (val_loss={val_loss:.4f})")
            else:
                # No validation, save after each epoch
                torch.save(model.state_dict(), best_model_path)
            
            if (epoch + 1) % 10 == 0:
                logger.info(f"Epoch {epoch+1}/{config['epochs']}: train_loss={train_loss:.4f}")
        
        # Load best model
        model.load_state_dict(torch.load(best_model_path, map_location=device))
        model.eval()
        
        trained_models[prop] = {
            "model": model,
            "config": config,
            "path": str(best_model_path),
            "best_val_loss": best_val_loss,
        }
        
        logger.info(f"Training complete for {prop}. Best val_loss: {best_val_loss:.4f}")
    
    return trained_models


def _collate_fn(batch):
    """Collate function for DataLoader."""
    data_list, targets = zip(*batch)
    batch_data = Batch.from_data_list(data_list)
    batch_targets = torch.stack(targets)
    return batch_data, batch_targets


def _create_mock_models(model_configs: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """Create mock models when PyTorch is not available."""
    models_dir = Path("data/models")
    models_dir.mkdir(parents=True, exist_ok=True)
    
    mock_models = {}
    for prop, config in model_configs.items():
        mock_models[prop] = {
            "model": None,
            "config": config,
            "path": str(models_dir / f"{prop}_best.pt"),
            "best_val_loss": 0.5,  # Mock
        }
    
    return mock_models


if __name__ == "__main__":
    # This would be called from train_ml_models.py
    pass

