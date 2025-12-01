"""
Model Registry - Manage ML models

Handles model registration, loading, and metadata.
"""

import json
import os
from typing import Dict, List, Optional, Any
from pathlib import Path
import logging

# Optional PyTorch imports - lazy import
TORCH_AVAILABLE = None
torch = None

def _ensure_torch():
    """Lazy import of PyTorch."""
    global TORCH_AVAILABLE, torch
    if TORCH_AVAILABLE is not None:
        return TORCH_AVAILABLE
    
    try:
        import torch as _torch
        torch = _torch
        TORCH_AVAILABLE = True
        return True
    except (ImportError, OSError) as e:
        TORCH_AVAILABLE = False
        torch = None
        logging.warning(f"PyTorch not available: {e}. Model registry will use mock mode.")
        return False

try:
    from .gat_model import AttentionGNN, create_model
except (ImportError, OSError):
    AttentionGNN = None
    create_model = None


class ModelRegistry:
    """
    Registry for managing ML models.
    """
    
    def __init__(self, registry_path: str = "models/registry.json"):
        self.registry_path = registry_path
        self.models: Dict[str, Dict[str, Any]] = {}
        self.loaded_models: Dict[str, Any] = {}
        _ensure_torch()  # Try to load PyTorch
        if TORCH_AVAILABLE and torch is not None:
            self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        else:
            self.device = None
        
        # Load registry if exists
        if os.path.exists(registry_path):
            self.load_registry()
        else:
            # Initialize with defaults
            self._init_defaults()
    
    def _init_defaults(self):
        """Initialize with default models."""
        self.register(
            model_id="attention-gnn-default",
            model_type="attention-gnn",
            path="models/attention_gnn_default.pt",
            metadata={
                "node_feat_dim": 64,
                "edge_feat_dim": 8,
                "hidden_dim": 128,
                "out_dim": 1,
                "num_layers": 3,
                "heads": 4,
            },
            description="Default Attention GNN for property prediction",
        )
    
    def register(
        self,
        model_id: str,
        model_type: str,
        path: str,
        metadata: Dict[str, Any],
        description: Optional[str] = None,
        is_default: bool = False,
    ):
        """
        Register a model.
        
        Args:
            model_id: Unique model identifier
            model_type: Type of model ("attention-gnn", "gnn", "classical", etc.)
            path: Path to model checkpoint
            metadata: Model configuration/metadata
            description: Human-readable description
            is_default: Whether this is the default model
        """
        self.models[model_id] = {
            "id": model_id,
            "type": model_type,
            "path": path,
            "metadata": metadata,
            "description": description,
            "is_default": is_default,
        }
        
        # If this is default, unset others
        if is_default:
            for mid in self.models:
                if mid != model_id:
                    self.models[mid]["is_default"] = False
        
        self.save_registry()
    
    def get_model(self, model_id: str) -> Optional[Dict[str, Any]]:
        """Get model metadata."""
        return self.models.get(model_id)
    
    def list_models(self) -> List[Dict[str, Any]]:
        """List all registered models."""
        return list(self.models.values())
    
    def get_default_model_id(self) -> Optional[str]:
        """Get default model ID."""
        for model_id, model_info in self.models.items():
            if model_info.get("is_default"):
                return model_id
        
        # Return first model if no default
        if self.models:
            return list(self.models.keys())[0]
        
        return None
    
    def load_model(self, model_id: str) -> Optional[Any]:
        """
        Load model from checkpoint.
        
        Returns:
            Loaded model instance or None if failed
        """
        # Check if already loaded
        if model_id in self.loaded_models:
            return self.loaded_models[model_id]
        
        model_info = self.get_model(model_id)
        if not model_info:
            return None
        
        model_type = model_info["type"]
        path = model_info["path"]
        metadata = model_info["metadata"]
        
        if not os.path.exists(path):
            print(f"Warning: Model checkpoint not found at {path}")
            return None
        
        try:
            if model_type == "attention-gnn":
                # Create model with metadata
                model = create_model(
                    node_feat_dim=metadata.get("node_feat_dim", 64),
                    edge_feat_dim=metadata.get("edge_feat_dim", 8),
                    hidden_dim=metadata.get("hidden_dim", 128),
                    out_dim=metadata.get("out_dim", 1),
                    num_layers=metadata.get("num_layers", 3),
                    heads=metadata.get("heads", 4),
                )
                
                # Load weights
                model.load_state_dict(torch.load(path, map_location=self.device))
                model.eval()
                model.to(self.device)
                
                self.loaded_models[model_id] = model
                return model
            
            elif model_type == "gnn":
                # Similar for standard GNN
                # TODO: Implement
                return None
            
            else:
                print(f"Unknown model type: {model_type}")
                return None
        
        except Exception as e:
            print(f"Failed to load model {model_id}: {e}")
            return None
    
    def unload_model(self, model_id: str):
        """Unload model from memory."""
        if model_id in self.loaded_models:
            del self.loaded_models[model_id]
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
    
    def save_registry(self):
        """Save registry to file."""
        os.makedirs(os.path.dirname(self.registry_path), exist_ok=True)
        with open(self.registry_path, 'w') as f:
            json.dump(self.models, f, indent=2)
    
    def load_registry(self):
        """Load registry from file."""
        try:
            with open(self.registry_path, 'r') as f:
                self.models = json.load(f)
        except Exception as e:
            print(f"Failed to load registry: {e}")
            self.models = {}
            self._init_defaults()
    
    def remove_model(self, model_id: str) -> bool:
        """Remove model from registry."""
        if model_id in self.models:
            self.unload_model(model_id)
            del self.models[model_id]
            self.save_registry()
            return True
        return False

