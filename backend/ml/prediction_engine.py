"""
Prediction Engine - Unified API for molecular property prediction

Supports multiple model types (classical, GNN, Attention-GNN) with fallbacks.
"""

from typing import Dict, List, Optional, Tuple, Any
import torch
import numpy as np
from torch_geometric.data import Data, Batch
import logging

from .featurize import featurize_smiles, featurize_json, canonicalize_smiles
from .registry import ModelRegistry
from .gat_model import AttentionGNN

logger = logging.getLogger(__name__)


class PredictionResult:
    """Structured prediction result."""
    
    def __init__(
        self,
        predictions: Dict[str, float],
        confidence: Optional[Dict[str, float]] = None,
        attention_weights: Optional[Dict[str, List[float]]] = None,
        edge_index: Optional[List[List[int]]] = None,
        node_mapping: Optional[Dict[int, str]] = None,
        model_id: Optional[str] = None,
        warnings: Optional[List[str]] = None,
    ):
        self.predictions = predictions
        self.confidence = confidence or {}
        self.attention_weights = attention_weights
        self.edge_index = edge_index
        self.node_mapping = node_mapping
        self.model_id = model_id
        self.warnings = warnings or []
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        result = {
            "predictions": self.predictions,
            "model_id": self.model_id,
        }
        
        if self.confidence:
            result["confidence"] = self.confidence
        
        if self.attention_weights:
            result["attentions"] = self.attention_weights
            if self.edge_index:
                result["edge_index"] = self.edge_index
            if self.node_mapping:
                result["node_mapping"] = self.node_mapping
        
        if self.warnings:
            result["warnings"] = self.warnings
        
        return result


class PredictionEngine:
    """
    Unified prediction engine supporting multiple model types.
    """
    
    def __init__(self, model_registry: Optional[ModelRegistry] = None):
        self.registry = model_registry or ModelRegistry()
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.loaded_models: Dict[str, Any] = {}
        self.feature_cache: Dict[str, Data] = {}
    
    def predict(
        self,
        input_data: Dict[str, Any],
        model_id: Optional[str] = None,
        properties: Optional[List[str]] = None,
        return_attention: bool = False,
        attention_layer: str = "last",
    ) -> PredictionResult:
        """
        Main prediction method.
        
        Args:
            input_data: Can be {"smiles": "..."} or {"graph": {...}} or {"molecule": {...}}
            model_id: Model to use (default: best available)
            properties: List of properties to predict (default: all)
            return_attention: Whether to return attention weights
            attention_layer: "last", "all", or layer index
        
        Returns:
            PredictionResult with predictions and optional attention
        """
        # Determine input type and featurize
        if "smiles" in input_data:
            smiles = input_data["smiles"]
            atom_order = input_data.get("node_order")
            
            # Check cache
            cache_key = f"smiles:{canonicalize_smiles(smiles)}"
            if cache_key in self.feature_cache:
                data = self.feature_cache[cache_key]
                node_mapping = {}  # Would need to store this too
            else:
                data, node_mapping = featurize_smiles(smiles, atom_order)
                self.feature_cache[cache_key] = data
            
        elif "graph" in input_data:
            data, node_mapping = featurize_json(input_data["graph"])
        elif "molecule" in input_data:
            data, node_mapping = featurize_json(input_data["molecule"])
        else:
            raise ValueError("Invalid input format. Expected 'smiles', 'graph', or 'molecule'")
        
        # Select model
        if not model_id:
            model_id = self.registry.get_default_model_id()
        
        if not model_id:
            raise ValueError("No model available")
        
        # Load model if not already loaded
        if model_id not in self.loaded_models:
            model = self.registry.load_model(model_id)
            if model:
                model.eval()
                model.to(self.device)
                self.loaded_models[model_id] = model
            else:
                raise ValueError(f"Failed to load model: {model_id}")
        
        model = self.loaded_models[model_id]
        model_info = self.registry.get_model(model_id)
        
        # Run inference
        batch = Batch.from_data_list([data]).to(self.device)
        
        with torch.no_grad():
            if isinstance(model, AttentionGNN) and return_attention:
                # Attention model
                preds, attentions = model(
                    batch.x,
                    batch.edge_index,
                    edge_attr=getattr(batch, 'edge_attr', None),
                    batch=batch.batch,
                    return_attentions=True,
                )
                
                # Process attentions
                attention_dict = self._process_attentions(
                    attentions,
                    attention_layer,
                    data.edge_index,
                )
                
                # Compute node importance
                node_importance = self._compute_node_importance(
                    attention_dict.get("layer_last", []),
                    data.edge_index,
                )
                
            else:
                # Standard model
                preds = model(
                    batch.x,
                    batch.edge_index,
                    edge_attr=getattr(batch, 'edge_attr', None),
                    batch=batch.batch,
                )
                attention_dict = None
                node_importance = None
        
        # Convert predictions to dict
        pred_values = preds.cpu().numpy().flatten()
        
        # Map to property names
        if properties:
            predictions = {prop: float(pred_values[i]) for i, prop in enumerate(properties) if i < len(pred_values)}
        else:
            # Default property names
            default_props = ["logP", "toxicity", "solubility", "molecularWeight"]
            predictions = {prop: float(pred_values[i]) if i < len(pred_values) else 0.0 
                          for i, prop in enumerate(default_props)}
        
        # Build result
        result = PredictionResult(
            predictions=predictions,
            model_id=model_id,
        )
        
        if attention_dict:
            result.attention_weights = attention_dict
            result.edge_index = data.edge_index.cpu().numpy().tolist()
            result.node_mapping = node_mapping
            if node_importance:
                result.attention_weights["node_importance"] = node_importance
        
        # Sanity checks
        warnings = self._sanity_check(predictions, data)
        result.warnings = warnings
        
        return result
    
    def _process_attentions(
        self,
        attentions: List[torch.Tensor],
        layer_selection: str,
        edge_index: torch.Tensor,
    ) -> Dict[str, List[float]]:
        """Process and normalize attention weights."""
        from .attention_utils import minmax_normalize
        
        result = {}
        
        if layer_selection == "all":
            for i, att in enumerate(attentions):
                normalized = minmax_normalize(att.cpu().numpy().tolist())
                result[f"layer_{i}"] = normalized
        elif layer_selection == "last":
            if attentions:
                last_att = attentions[-1]
                normalized = minmax_normalize(last_att.cpu().numpy().tolist())
                result["layer_last"] = normalized
        else:
            # Layer index
            try:
                idx = int(layer_selection)
                if 0 <= idx < len(attentions):
                    att = attentions[idx]
                    normalized = minmax_normalize(att.cpu().numpy().tolist())
                    result[f"layer_{idx}"] = normalized
            except ValueError:
                pass
        
        return result
    
    def _compute_node_importance(
        self,
        edge_attentions: List[float],
        edge_index: torch.Tensor,
    ) -> List[float]:
        """Compute node importance by aggregating edge attentions."""
        if not edge_attentions or edge_index.size(1) == 0:
            return []
        
        num_nodes = edge_index.max().item() + 1
        node_scores = np.zeros(num_nodes)
        
        edge_idx_np = edge_index.cpu().numpy()
        
        for e in range(len(edge_attentions)):
            if e < edge_idx_np.shape[1]:
                u = int(edge_idx_np[0, e])
                v = int(edge_idx_np[1, e])
                att = edge_attentions[e]
                node_scores[u] += att
                node_scores[v] += att
        
        # Normalize
        max_score = node_scores.max()
        if max_score > 0:
            node_scores = node_scores / max_score
        
        return node_scores.tolist()
    
    def _sanity_check(
        self,
        predictions: Dict[str, float],
        data: Data,
    ) -> List[str]:
        """Perform sanity checks on predictions."""
        warnings = []
        
        # Check for extreme values
        for prop, value in predictions.items():
            if prop == "logP":
                if value < -10 or value > 20:
                    warnings.append(f"logP value {value:.2f} is outside typical range [-10, 20]")
            elif prop == "toxicity":
                if value < 0 or value > 1:
                    warnings.append(f"Toxicity value {value:.2f} should be in [0, 1]")
            elif prop == "molecularWeight":
                if value < 0 or value > 10000:
                    warnings.append(f"Molecular weight {value:.2f} seems unrealistic")
        
        # Check molecule size
        if data.num_nodes > 500:
            warnings.append(f"Large molecule ({data.num_nodes} atoms) - predictions may be less reliable")
        
        return warnings
    
    def predict_batch(
        self,
        inputs: List[Dict[str, Any]],
        model_id: Optional[str] = None,
        batch_size: int = 32,
    ) -> List[PredictionResult]:
        """Batch prediction for multiple molecules."""
        results = []
        
        for i in range(0, len(inputs), batch_size):
            batch_inputs = inputs[i:i + batch_size]
            
            # Featurize batch
            data_list = []
            for inp in batch_inputs:
                if "smiles" in inp:
                    data, _ = featurize_smiles(inp["smiles"], inp.get("node_order"))
                elif "graph" in inp:
                    data, _ = featurize_json(inp["graph"])
                else:
                    continue
                data_list.append(data)
            
            if not data_list:
                continue
            
            # Batch inference
            batch = Batch.from_data_list(data_list).to(self.device)
            
            # Load model
            if not model_id:
                model_id = self.registry.get_default_model_id()
            
            if model_id not in self.loaded_models:
                model = self.registry.load_model(model_id)
                if model:
                    model.eval()
                    model.to(self.device)
                    self.loaded_models[model_id] = model
            
            model = self.loaded_models[model_id]
            
            with torch.no_grad():
                preds = model(
                    batch.x,
                    batch.edge_index,
                    edge_attr=getattr(batch, 'edge_attr', None),
                    batch=batch.batch,
                )
            
            # Split predictions
            preds_np = preds.cpu().numpy()
            num_per_graph = [data.num_nodes for data in data_list]
            # This is simplified - actual splitting requires proper batching logic
            
            for j, inp in enumerate(batch_inputs):
                if j < len(preds_np):
                    predictions = {"logP": float(preds_np[j][0])}
                    results.append(PredictionResult(
                        predictions=predictions,
                        model_id=model_id,
                    ))
        
        return results
    
    def clear_cache(self):
        """Clear feature cache."""
        self.feature_cache.clear()
    
    def unload_model(self, model_id: str):
        """Unload model from memory."""
        if model_id in self.loaded_models:
            del self.loaded_models[model_id]
            if torch.cuda.is_available():
                torch.cuda.empty_cache()

