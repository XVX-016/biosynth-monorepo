"""
Training script for PropertyPredictor model
"""
import os
import sys
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from backend.ai.property_predictor import create_model
from backend.ai.featurizer import featurize_smiles


class MoleculeDataset(Dataset):
    """Dataset for molecular property prediction"""
    
    def __init__(self, csv_path: str):
        """
        Initialize dataset from CSV file
        
        Expected CSV format:
        smiles,stability,toxicity,solubility,bioavailability,novelty
        """
        self.df = pd.read_csv(csv_path)
        self.smiles = self.df['smiles'].values
        self.properties = self.df[['stability', 'toxicity', 'solubility', 'bioavailability', 'novelty']].values
    
    def __len__(self):
        return len(self.df)
    
    def __getitem__(self, idx):
        smiles = self.smiles[idx]
        properties = self.properties[idx]
        
        # Featurize SMILES
        features = featurize_smiles(smiles)
        if features is None:
            # Invalid SMILES - return zeros
            features = [0.0] * 2048
        
        return {
            'features': torch.tensor(features, dtype=torch.float32),
            'properties': torch.tensor(properties, dtype=torch.float32)
        }


def generate_dummy_dataset(output_path: str, num_samples: int = 50):
    """
    Generate a dummy dataset for training
    
    Args:
        output_path: Path to save CSV file
        num_samples: Number of samples to generate
    """
    # Common SMILES patterns for dummy data
    smiles_list = [
        'C', 'CC', 'CCC', 'CCO', 'CCN', 'CCOC', 'CC(=O)O', 'CCN(CC)CC',
        'c1ccccc1', 'CCc1ccccc1', 'CC(C)C', 'CC(C)(C)C', 'CCOC(=O)C',
        'CCN', 'CCCC', 'CCCCC', 'CCCCCC', 'CCCCCCC', 'CCCCCCCC',
        'C1CCCCC1', 'C1CCCC1', 'C1=CC=CC=C1', 'C1=CC=CC=C1C',
        'CC(C)CC', 'CC(C)(C)CC', 'CCOCC', 'CCNCC', 'CC(=O)N',
        'CC(=O)OC', 'CC(=O)NC', 'CCN(C)C', 'CCOCCO', 'CCNCCN',
        'C1CCC1', 'C1CCCC1', 'C1CCCCC1', 'C1CCCCCC1', 'C1CCCCCCC1',
        'CC(C)(C)C', 'CC(C)(C)CC', 'CC(C)(C)CCC', 'CC(C)(C)CCCC',
        'CCO', 'CCOC', 'CCOCC', 'CCOCCC', 'CCOCCCC',
        'CCN', 'CCNC', 'CCNCC', 'CCNCCC', 'CCNCCCC',
        'CC(=O)O', 'CC(=O)OC', 'CC(=O)OCC', 'CC(=O)OCCC'
    ]
    
    # Generate data
    data = []
    for i in range(num_samples):
        smiles = smiles_list[i % len(smiles_list)]
        
        # Generate dummy properties (random but realistic ranges)
        stability = np.random.uniform(0.3, 0.9)
        toxicity = np.random.uniform(0.1, 0.7)
        solubility = np.random.uniform(0.2, 0.8)
        bioavailability = np.random.uniform(0.2, 0.9)
        novelty = np.random.uniform(0.1, 0.6)
        
        data.append({
            'smiles': smiles,
            'stability': stability,
            'toxicity': toxicity,
            'solubility': solubility,
            'bioavailability': bioavailability,
            'novelty': novelty
        })
    
    # Create DataFrame and save
    df = pd.DataFrame(data)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df.to_csv(output_path, index=False)
    print(f"Generated dummy dataset with {num_samples} samples at {output_path}")


def train_model(
    csv_path: str = 'backend/data/molecules.csv',
    epochs: int = 50,
    batch_size: int = 8,
    learning_rate: float = 0.001,
    output_path: str = 'backend/weights/property_predictor.pt'
):
    """
    Train PropertyPredictor model
    
    Args:
        csv_path: Path to training CSV file
        epochs: Number of training epochs
        batch_size: Batch size for training
        learning_rate: Learning rate for optimizer
        output_path: Path to save trained model
    """
    # Generate dummy dataset if it doesn't exist
    if not os.path.exists(csv_path):
        print(f"Dataset not found at {csv_path}. Generating dummy dataset...")
        generate_dummy_dataset(csv_path)
    
    # Create dataset and dataloader
    dataset = MoleculeDataset(csv_path)
    dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)
    
    # Create model
    model = create_model()
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)
    
    # Loss and optimizer
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)
    
    # Training loop
    model.train()
    for epoch in range(epochs):
        total_loss = 0.0
        num_batches = 0
        
        progress_bar = tqdm(dataloader, desc=f'Epoch {epoch+1}/{epochs}')
        for batch in progress_bar:
            features = batch['features'].to(device)
            properties = batch['properties'].to(device)
            
            # Forward pass
            optimizer.zero_grad()
            predictions = model(features)
            loss = criterion(predictions, properties)
            
            # Backward pass
            loss.backward()
            optimizer.step()
            
            total_loss += loss.item()
            num_batches += 1
            
            # Update progress bar
            progress_bar.set_postfix({'loss': f'{loss.item():.4f}'})
        
        avg_loss = total_loss / num_batches
        print(f'Epoch {epoch+1}/{epochs} - Average Loss: {avg_loss:.4f}')
    
    # Save model
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    torch.save(model.state_dict(), output_path)
    print(f"Model saved to {output_path}")
    
    return model


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Train PropertyPredictor model')
    parser.add_argument('--csv', type=str, default='backend/data/molecules.csv',
                       help='Path to training CSV file')
    parser.add_argument('--epochs', type=int, default=50,
                       help='Number of training epochs')
    parser.add_argument('--batch-size', type=int, default=8,
                       help='Batch size for training')
    parser.add_argument('--lr', type=float, default=0.001,
                       help='Learning rate')
    parser.add_argument('--output', type=str, default='backend/weights/property_predictor.pt',
                       help='Path to save trained model')
    
    args = parser.parse_args()
    
    train_model(
        csv_path=args.csv,
        epochs=args.epochs,
        batch_size=args.batch_size,
        learning_rate=args.lr,
        output_path=args.output
    )

