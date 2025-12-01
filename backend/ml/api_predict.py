"""
FastAPI endpoint for attention-based predictions

Place this in your backend/ml/ directory and integrate with your FastAPI app.
"""

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
import torch
from typing import List, Optional
from gat_model import AttentionGNN, create_model
from torch_geometric.data import Data, Batch

app = FastAPI()

# Load model (adjust path and dimensions as needed)
MODEL_PATH = "checkpoints/best_model.pt"
device = 'cuda' if torch.cuda.is_available() else 'cpu'

# Initialize model (adjust dimensions to match your trained model)
node_feat_dim = 64
edge_feat_dim = 8
model = create_model(
    node_feat_dim=node_feat_dim,
    edge_feat_dim=edge_feat_dim,
    hidden_dim=128,
    out_dim=1,
    num_layers=3,
    heads=4,
)

# Load weights
try:
    model.load_state_dict(torch.load(MODEL_PATH, map_location=device))
    model.eval()
    model.to(device)
    print(f"Model loaded from {MODEL_PATH}")
except FileNotFoundError:
    print(f"Warning: Model not found at {MODEL_PATH}, using untrained model")


class GraphFeature(BaseModel):
    """Graph feature representation from frontend."""
    nodes: List[List[float]]  # [N, node_feat_dim]
    edges: List[List[int]]  # [E, 2]
    edgeFeatures: Optional[List[List[float]]] = None  # [E, edge_feat_dim]


class PredictRequest(BaseModel):
    graph: GraphFeature
    modelId: Optional[str] = "attention-gnn"
    returnAttention: bool = True
    layerIndex: Optional[int] = None  # Which layer to return (None = last)


def featurize_graph(graph: GraphFeature) -> Data:
    """
    Convert frontend graph features to PyTorch Geometric Data object.
    """
    # Node features
    x = torch.tensor(graph.nodes, dtype=torch.float)
    
    # Edge indices
    if len(graph.edges) > 0:
        edge_index = torch.tensor(graph.edges, dtype=torch.long).t().contiguous()
    else:
        edge_index = torch.empty((2, 0), dtype=torch.long)
    
    # Edge features
    edge_attr = None
    if graph.edgeFeatures and len(graph.edgeFeatures) > 0:
        edge_attr = torch.tensor(graph.edgeFeatures, dtype=torch.float)
    
    return Data(x=x, edge_index=edge_index, edge_attr=edge_attr)


@app.post("/api/ml/predict")
async def predict(req: PredictRequest):
    """
    Predict molecular properties with attention weights.
    """
    try:
        # Convert to PyG Data
        data = featurize_graph(req.graph)
        batch = Batch.from_data_list([data]).to(device)
        
        # Forward pass
        with torch.no_grad():
            preds, attentions = model(
                batch.x,
                batch.edge_index,
                edge_attr=getattr(batch, 'edge_attr', None),
                batch=batch.batch,
                return_attentions=True,
            )
        
        # Extract predictions
        prediction = preds.cpu().numpy().tolist()[0]
        
        # Select layer (default: last)
        layer_idx = req.layerIndex if req.layerIndex is not None else len(attentions) - 1
        if layer_idx >= len(attentions):
            layer_idx = len(attentions) - 1
        
        attention_layer = attentions[layer_idx] if attentions else None
        
        # Compute node importance (sum of incident edge attentions)
        node_importance = None
        if attention_layer is not None:
            num_nodes = data.x.size(0)
            node_importance = torch.zeros(num_nodes)
            edge_index = data.edge_index.cpu()
            
            for e in range(edge_index.size(1)):
                u = edge_index[0, e].item()
                v = edge_index[1, e].item()
                att = attention_layer[e].item()
                node_importance[u] += att
                node_importance[v] += att
            
            # Normalize
            max_imp = node_importance.max()
            if max_imp > 0:
                node_importance = node_importance / max_imp
        
        # Return results
        result = {
            "prediction": prediction,
            "edge_index": batch.edge_index.cpu().numpy().tolist(),
        }
        
        if req.returnAttention and attention_layer is not None:
            result["attentions"] = [attention_layer.numpy().tolist()]
            if node_importance is not None:
                result["node_importance"] = node_importance.numpy().tolist()
        
        return result
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/ml/predict/batch")
async def predict_batch(features: List[GraphFeature], modelId: Optional[str] = "attention-gnn"):
    """
    Batch prediction endpoint.
    """
    try:
        # Convert all graphs
        data_list = [featurize_graph(g) for g in features]
        batch = Batch.from_data_list(data_list).to(device)
        
        # Forward pass
        with torch.no_grad():
            preds, attentions = model(
                batch.x,
                batch.edge_index,
                edge_attr=getattr(batch, 'edge_attr', None),
                batch=batch.batch,
                return_attentions=True,
            )
        
        # Split predictions by graph
        predictions = []
        for i in range(len(features)):
            pred = preds[i].cpu().numpy().tolist()
            # Extract attention for this graph (simplified - would need proper splitting)
            predictions.append({
                "prediction": pred,
                "attentions": attentions[-1].numpy().tolist() if attentions else [],
            })
        
        return {"predictions": predictions}
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

