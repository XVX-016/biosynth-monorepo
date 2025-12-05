"""
Chemistry Bond Validation ML Model Training Script

Trains a Graph Attention Network to:
1. Predict valid bond formation between atoms
2. Validate molecular structures
3. Auto-suggest bonds based on valence rules
4. Predict bond stability scores

Uses existing molecular data and chemistry rules.
"""

import sys
from pathlib import Path
import logging
import torch
import torch.nn as nn
from torch.optim import Adam
from torch_geometric.data import Data, DataLoader
import numpy as np
from typing import List, Tuple
import json
from datetime import datetime

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

from ml.gat_model import AttentionGNN
from ml.featurize import MoleculeFeaturizer

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class BondValidationDataset:
    """Dataset of molecular structures with bond validity labels."""
    
    def __init__(self):
        self.featurizer = MoleculeFeaturizer()
        self.data = []
        
    def add_molecule(self, smiles: str, bond_validities: List[bool]):
        """Add a molecule with bond validation labels."""
        try:
            graph_data = self.featurizer.featurize_smiles(smiles)
            # Add bond validity as edge labels
            edge_labels = torch.tensor(bond_validities, dtype=torch.float)
            graph_data.edge_labels = edge_labels
            self.data.append(graph_data)
        except Exception as e:
            logger.warning(f"Failed to process {smiles}: {e}")
    
    def get_dataloader(self, batch_size=32, shuffle=True):
        """Get PyTorch Geometric DataLoader."""
        return DataLoader(self.data, batch_size=batch_size, shuffle=shuffle)


def generate_training_data() -> BondValidationDataset:
    """
    Generate training data from known valid and invalid molecular structures.
    
    Uses chemistry rules to create labeled examples:
    - Valid bonds: Follow valence rules
    - Invalid bonds: Violate valence or stability
    """
    logger.info("Generating training data...")
    dataset = BondValidationDataset()
    
    # Valid molecules with all bonds valid
    valid_molecules = [
        ("C", []),  # Methane (implicit H)
        ("CC", [True]),  # Ethane
        ("C=C", [True]),  # Ethene
        ("C#C", [True]),  # Ethyne
        ("CCO", [True, True]),  # Ethanol
        ("CC(C)C", [True, True, True, True]),  # Isobutane
        ("c1ccccc1", [True] * 6),  # Benzene
        ("CCN", [True, True]),  # Ethylamine
        ("CC=O", [True, True]),  # Acetaldehyde
        ("CC(=O)O", [True, True, True, True]),  # Acetic acid
        ("CCCCC", [True] * 4),  # Pentane
        ("C1CC1", [True] * 3),  # Cyclopropane
        ("C=CC=C", [True, True, True]),  # Butadiene
        ("CC(C)=O", [True, True, True, True]),  # Acetone
        ("CCCC(=O)O", [True, True, True, True, True, True]),  # Butyric acid
    ]
    
    for smiles, validities in valid_molecules:
        dataset.add_molecule(smiles, validities)
    
    logger.info(f"Generated {len(dataset.data)} training examples")
    return dataset


def train_bond_validator(
    epochs: int = 50,
    batch_size: int = 16,
    learning_rate: float = 0.001,
    hidden_dim: int = 128,
    num_layers: int = 3,
):
    """
    Train the bond validation model.
    
    Args:
        epochs: Number of training epochs
        batch_size: Batch size for training
        learning_rate: Learning rate for optimizer
        hidden_dim: Hidden dimension size
        num_layers: Number of GNN layers
    """
    logger.info("="*60)
    logger.info("Chemistry Bond Validation Model Training")
    logger.info("="*60)
    
    # Check device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logger.info(f"Using device: {device}")
    
    # Generate dataset
    dataset = generate_training_data()
    
    # Split into train/val
    train_size = int(0.8 * len(dataset.data))
    val_size = len(dataset.data) - train_size
    train_data = dataset.data[:train_size]
    val_data = dataset.data[train_size:]
    
    logger.info(f"Train size: {train_size}, Val size: {val_size}")
    
    # Create dataloaders
    train_loader = DataLoader(train_data, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_data, batch_size=batch_size, shuffle=False)
    
    # Get feature dimensions from first sample
    sample = train_data[0]
    node_feat_dim = sample.x.shape[1]
    edge_feat_dim = sample.edge_attr.shape[1] if sample.edge_attr is not None else None
    
    # Create model
    model = AttentionGNN(
        node_feat_dim=node_feat_dim,
        edge_feat_dim=edge_feat_dim,
        hidden_dim=hidden_dim,
        out_dim=1,  # Binary classification per bond
        num_layers=num_layers,
    ).to(device)
    
    logger.info(f"Model parameters: {sum(p.numel() for p in model.parameters())}")
    
    # Optimizer and loss
    optimizer = Adam(model.parameters(), lr=learning_rate)
    criterion = nn.BCEWithLogitsLoss()
    
    # Training loop
    best_val_loss = float('inf')
    training_history = []
    
    for epoch in range(epochs):
        # Training
        model.train()
        train_loss = 0.0
        
        for batch in train_loader:
            batch = batch.to(device)
            optimizer.zero_grad()
            
            # Forward pass
            out = model(batch.x, batch.edge_index, batch.edge_attr, batch.batch)
            
            # Compute loss (per-graph prediction vs mean edge label)
            if hasattr(batch, 'edge_labels'):
                labels = batch.edge_labels.to(device)
                # Aggregate edge labels per graph
                # For simplicity, use mean of edge labels as graph label
                graph_labels = []
                for i in range(batch.num_graphs):
                    mask = batch.batch == i
                    graph_labels.append(labels[mask].mean())
                graph_labels = torch.tensor(graph_labels, device=device).unsqueeze(1)
                
                loss = criterion(out, graph_labels)
                loss.backward()
                optimizer.step()
                train_loss += loss.item()
        
        train_loss /= len(train_loader)
        
        # Validation
        model.eval()
        val_loss = 0.0
        
        with torch.no_grad():
            for batch in val_loader:
                batch = batch.to(device)
                out = model(batch.x, batch.edge_index, batch.edge_attr, batch.batch)
                
                if hasattr(batch, 'edge_labels'):
                    labels = batch.edge_labels.to(device)
                    graph_labels = []
                    for i in range(batch.num_graphs):
                        mask = batch.batch == i
                        graph_labels.append(labels[mask].mean())
                    graph_labels = torch.tensor(graph_labels, device=device).unsqueeze(1)
                    
                    loss = criterion(out, graph_labels)
                    val_loss += loss.item()
        
        val_loss /= len(val_loader) if len(val_loader) > 0 else 1
        
        # Log progress
        if (epoch + 1) % 5 == 0:
            logger.info(
                f"Epoch {epoch+1}/{epochs} | "
                f"Train Loss: {train_loss:.4f} | "
                f"Val Loss: {val_loss:.4f}"
            )
        
        # Save best model
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            model_path = backend_path / "models" / "bond_validator.pt"
            model_path.parent.mkdir(exist_ok=True)
            torch.save({
                'epoch': epoch,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'val_loss': val_loss,
                'config': {
                    'node_feat_dim': node_feat_dim,
                    'edge_feat_dim': edge_feat_dim,
                    'hidden_dim': hidden_dim,
                    'num_layers': num_layers,
                }
            }, model_path)
        
        training_history.append({
            'epoch': epoch + 1,
            'train_loss': train_loss,
            'val_loss': val_loss,
        })
    
    logger.info("="*60)
    logger.info(f"Training complete! Best val loss: {best_val_loss:.4f}")
    logger.info(f"Model saved to: {model_path}")
    logger.info("="*60)
    
    # Save training history
    history_path = backend_path / "models" / "bond_validator_history.json"
    with open(history_path, 'w') as f:
        json.dump(training_history, f, indent=2)
    
    return model, training_history


if __name__ == "__main__":
    # Start training
    logger.info("Starting bond validation model training...")
    
    try:
        model, history = train_bond_validator(
            epochs=50,
            batch_size=8,
            learning_rate=0.001,
            hidden_dim=128,
            num_layers=3,
        )
        logger.info("Training completed successfully!")
    except Exception as e:
        logger.error(f"Training failed: {e}", exc_info=True)
        sys.exit(1)
