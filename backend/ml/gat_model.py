"""
AttentionGNN - Graph Attention Network model for molecular property prediction

This is a backend Python implementation using PyTorch Geometric.
Place this in your backend/ml/ directory.
"""

import torch
import torch.nn.functional as F
from torch import nn
from torch_geometric.nn import GATConv, global_mean_pool


class AttentionGNN(nn.Module):
    """
    Graph Attention Network for molecular property prediction.
    
    Returns attention weights for explainability.
    """
    
    def __init__(
        self,
        node_feat_dim: int,
        edge_feat_dim: int = None,
        hidden_dim: int = 128,
        out_dim: int = 1,
        num_layers: int = 3,
        heads: int = 4,
        dropout: float = 0.1,
    ):
        super().__init__()
        
        self.node_feat_dim = node_feat_dim
        self.edge_feat_dim = edge_feat_dim
        self.hidden_dim = hidden_dim
        self.out_dim = out_dim
        self.num_layers = num_layers
        self.heads = heads
        self.dropout = dropout
        
        # Initial linear projection
        self.input_lin = nn.Linear(node_feat_dim, hidden_dim)
        
        # Stack GATConv layers
        self.convs = nn.ModuleList()
        for i in range(num_layers):
            in_dim = hidden_dim if i > 0 else hidden_dim
            # GATConv with edge features (requires PyG >= 2.0)
            # Last layer: use concat=False so output is hidden_dim // heads
            # But we'll add a projection to hidden_dim after
            self.convs.append(
                GATConv(
                    in_dim,
                    hidden_dim // heads,
                    heads=heads,
                    dropout=dropout,
                    edge_dim=edge_feat_dim,
                    concat=True if i < num_layers - 1 else False,
                )
            )
        
        # Projection layer: if last layer uses concat=False, output is hidden_dim // heads
        # We need to project it to hidden_dim for the MLP
        if num_layers > 0:
            # Last layer outputs hidden_dim // heads when concat=False
            # Or hidden_dim when concat=True
            # To be safe, we'll always project to hidden_dim
            self.final_proj = nn.Linear(hidden_dim // heads, hidden_dim)
        else:
            self.final_proj = nn.Identity()
        
        # Readout (global pooling)
        self.pool = global_mean_pool
        
        # Final MLP
        self.mlp = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, out_dim),
        )
    
    def forward(
        self,
        x: torch.Tensor,
        edge_index: torch.Tensor,
        edge_attr: torch.Tensor = None,
        batch: torch.Tensor = None,
        return_attentions: bool = False,
    ):
        """
        Forward pass.
        
        Args:
            x: Node features [N, node_feat_dim]
            edge_index: Edge indices [2, E]
            edge_attr: Edge features [E, edge_feat_dim] or None
            batch: Batch vector [N] or None
            return_attentions: If True, return attention weights
        
        Returns:
            predictions: [batch_size, out_dim] or [1, out_dim]
            attentions: List of [E] tensors (one per layer) if return_attentions=True
        """
        x = self.input_lin(x)
        attentions = []
        
        for conv in self.convs:
            # Get attention weights
            try:
                x, (edge_idx_out, alpha) = conv(
                    x, edge_index, edge_attr, return_attention_weights=True
                )
                
                # alpha shape: [E, heads] -> average across heads
                if alpha.dim() == 2:
                    alpha_edge = alpha.mean(dim=1)  # [E]
                else:
                    alpha_edge = alpha.view(-1)  # fallback
                
                attentions.append(alpha_edge.detach().cpu())
            except TypeError:
                # Older PyG versions - can't extract attention directly
                x = conv(x, edge_index, edge_attr)
                attentions.append(None)
            
            x = F.elu(x)
            x = F.dropout(x, p=self.dropout, training=self.training)
        
        # Global pooling
        if batch is not None:
            out = self.pool(x, batch)
        else:
            out = x.mean(dim=0, keepdim=True)
        
        # Project to hidden_dim if needed (when last layer used concat=False)
        # Check if output dimension matches expected hidden_dim
        if out.shape[-1] != self.hidden_dim:
            out = self.final_proj(out)
        
        # Final MLP
        out = self.mlp(out)
        
        if return_attentions:
            # Filter out None values
            attentions = [a for a in attentions if a is not None]
            return out, attentions
        
        return out


def create_model(
    node_feat_dim: int = 64,
    edge_feat_dim: int = 8,
    hidden_dim: int = 128,
    out_dim: int = 1,
    num_layers: int = 3,
    heads: int = 4,
) -> AttentionGNN:
    """Factory function to create model."""
    return AttentionGNN(
        node_feat_dim=node_feat_dim,
        edge_feat_dim=edge_feat_dim,
        hidden_dim=hidden_dim,
        out_dim=out_dim,
        num_layers=num_layers,
        heads=heads,
    )

