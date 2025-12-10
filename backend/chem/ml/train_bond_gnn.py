"""
Train a bond-prediction GNN. This is a skeleton showing data pipeline and model
using PyTorch Geometric. You must prepare a dataset of molecules where each
example contains node features and ground-truth bond adjacency / labels.

This script is intentionally minimal and annotated. Replace data loader with
RDKit-based featurizer that provides coordinates and atom features.
"""

import os
import torch
from torch.nn import BCEWithLogitsLoss
try:
    from torch_geometric.data import DataLoader
    from torch_geometric.nn import GCNConv, global_mean_pool
except ImportError:
    print("torch_geometric not installed. Skipping imports.")

class BondGNN(torch.nn.Module):
    def __init__(self, in_dim, hid=64):
        super().__init__()
        self.conv1 = GCNConv(in_dim, hid)
        self.conv2 = GCNConv(hid, hid)
        self.head = torch.nn.Linear(hid*2, 1)  # pairwise

    def forward(self, x, edge_index, batch, pairs):
        # x: [N, F], edge_index: [2, E]
        x = self.conv1(x, edge_index).relu()
        x = self.conv2(x, edge_index).relu()
        # compute pair embeddings
        emb_a = x[pairs[:,0]]
        emb_b = x[pairs[:,1]]
        out = self.head(torch.cat([emb_a, emb_b], dim=-1)).squeeze(-1)
        return out

# TODO: implement Dataset that yields (x, edge_index, batch, pairs, labels)
# For training loop example only:

def train():
    # placeholder dataset loader
    train_loader = []
    # in_dim placeholder
    model = BondGNN(in_dim=16)
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    loss_fn = BCEWithLogitsLoss()

    for epoch in range(1,21):
        model.train()
        total = 0.0
        for batch in train_loader:
            # batch.x, batch.edge_index, batch.batch, batch.pairs, batch.labels
            pred = model(batch.x, batch.edge_index, batch.batch, batch.pairs)
            loss = loss_fn(pred, batch.labels.float())
            opt.zero_grad(); loss.backward(); opt.step()
            total += loss.item()
        print('epoch', epoch, 'loss', total)

    # save a simple state dict
    os.makedirs('backend/chem/ml', exist_ok=True)
    torch.save(model.state_dict(), os.environ.get('BOND_MODEL_PATH','./backend/chem/ml/model.pth'))

if __name__ == '__main__':
    try:
        train()
    except Exception as e:
        print(f"Training failed (likely due to missing data/dependencies): {e}")
