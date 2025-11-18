"""
PyTorch neural network for molecular property prediction
"""
import torch
import torch.nn as nn
import torch.nn.functional as F


class PropertyPredictor(nn.Module):
    """
    Neural network for predicting molecular properties from Morgan fingerprints
    
    Architecture:
    - Input: 2048-bit Morgan fingerprint
    - Hidden layers: 1024 -> 512 -> 128
    - Output: 5 properties (stability, toxicity, solubility, bioavailability, novelty)
    """
    
    def __init__(self, input_dim: int = 2048, dropout: float = 0.3):
        super(PropertyPredictor, self).__init__()
        
        self.input_dim = input_dim
        self.output_dim = 5  # stability, toxicity, solubility, bioavailability, novelty
        
        # Hidden layers
        self.fc1 = nn.Linear(input_dim, 1024)
        self.fc2 = nn.Linear(1024, 512)
        self.fc3 = nn.Linear(512, 128)
        self.fc4 = nn.Linear(128, self.output_dim)
        
        self.dropout = nn.Dropout(dropout)
        
        # Initialize weights using Kaiming initialization
        self._initialize_weights()
    
    def _initialize_weights(self):
        """Initialize weights using Kaiming initialization"""
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.kaiming_normal_(m.weight, mode='fan_out', nonlinearity='relu')
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
    
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Forward pass
        
        Args:
            x: Input tensor of shape (batch_size, input_dim)
            
        Returns:
            Tensor of shape (batch_size, output_dim) with property predictions
        """
        # Input layer
        x = self.fc1(x)
        x = F.relu(x)
        x = self.dropout(x)
        
        # Hidden layer 1
        x = self.fc2(x)
        x = F.relu(x)
        x = self.dropout(x)
        
        # Hidden layer 2
        x = self.fc3(x)
        x = F.relu(x)
        x = self.dropout(x)
        
        # Output layer (no activation for regression)
        x = self.fc4(x)
        
        return x


def create_model(input_dim: int = 2048, dropout: float = 0.3) -> PropertyPredictor:
    """
    Factory function to create a PropertyPredictor model
    
    Args:
        input_dim: Dimension of input fingerprint (default: 2048)
        dropout: Dropout probability (default: 0.3)
        
    Returns:
        Initialized PropertyPredictor model
    """
    return PropertyPredictor(input_dim=input_dim, dropout=dropout)

