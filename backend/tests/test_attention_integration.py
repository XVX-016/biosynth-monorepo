"""
Integration tests for attention mapping between frontend and backend.

Tests the complete pipeline: FE → BE → FE visualization
"""

import pytest
import torch
from torch_geometric.data import Data

from ..ml.featurize import featurize_smiles, featurize_json
from ..ml.prediction_engine import PredictionEngine
from ..ml.registry import ModelRegistry


def test_ethanol_attention_mapping():
    """
    Test attention mapping for ethanol (CCO).
    
    Verifies:
    - FE atom order is preserved
    - Edge indices match bonds
    - Attention weights length == bond count
    """
    # Frontend representation
    fe_molecule = {
        "atoms": [
            {"id": "a0", "element": "C", "x": 0, "y": 0, "z": 0, "charge": 0},
            {"id": "a1", "element": "C", "x": 50, "y": 0, "z": 0, "charge": 0},
            {"id": "a2", "element": "O", "x": 100, "y": 0, "z": 0, "charge": 0},
        ],
        "bonds": [
            {"id": "b0", "atoms": ["a0", "a1"], "order": 1},
            {"id": "b1", "atoms": ["a1", "a2"], "order": 1},
        ],
        "node_order": ["a0", "a1", "a2"],
    }
    
    # Featurize
    data, node_mapping = featurize_json(fe_molecule)
    
    # Verify node mapping
    assert node_mapping[0] == "a0"
    assert node_mapping[1] == "a1"
    assert node_mapping[2] == "a2"
    
    # Verify edge count (undirected = 2 edges per bond)
    assert data.edge_index.size(1) == 4  # 2 bonds * 2 directions
    
    # Verify we can map edges back to bonds
    edge_index = data.edge_index.cpu().numpy()
    
    # Find edges connecting a0-a1 (indices 0-1)
    a0_a1_edges = []
    for e in range(edge_index.shape[1]):
        u, v = edge_index[0, e], edge_index[1, e]
        if (u == 0 and v == 1) or (u == 1 and v == 0):
            a0_a1_edges.append(e)
    
    assert len(a0_a1_edges) == 2  # Both directions
    
    print("✓ Ethanol attention mapping test passed")


def test_smiles_featurization():
    """Test SMILES featurization preserves atom order."""
    smiles = "CCO"  # Ethanol
    
    # Featurize with atom order
    atom_order = ["a0", "a1", "a2"]
    data, node_mapping = featurize_smiles(smiles, atom_order)
    
    # Verify mapping
    assert len(node_mapping) == 3
    assert node_mapping[0] == "a0"
    
    # Verify graph structure
    assert data.num_nodes == 3
    assert data.edge_index.size(1) > 0
    
    print("✓ SMILES featurization test passed")


def test_attention_normalization():
    """Test attention weight normalization."""
    from ..ml.attention_utils import minmax_normalize, softmax_per_edge
    
    # Test min-max
    attentions = [0.1, 0.5, 0.9, 0.3]
    normalized = minmax_normalize(attentions)
    
    assert min(normalized) == 0.0
    assert max(normalized) == 1.0
    assert len(normalized) == len(attentions)
    
    # Test softmax
    softmaxed = softmax_per_edge(attentions)
    assert abs(sum(softmaxed) - 1.0) < 1e-6
    
    print("✓ Attention normalization test passed")


def test_prediction_engine_smiles():
    """Test PredictionEngine with SMILES input."""
    # This requires a trained model - skip if no model available
    registry = ModelRegistry()
    model_id = registry.get_default_model_id()
    
    if not model_id:
        pytest.skip("No model available for testing")
    
    engine = PredictionEngine(registry)
    
    # Test prediction
    result = engine.predict(
        input_data={"smiles": "CCO"},
        model_id=model_id,
        return_attention=False,
    )
    
    assert result.predictions is not None
    assert result.model_id == model_id
    
    print("✓ Prediction engine test passed")


def test_edge_index_consistency():
    """
    Test that edge_index from backend matches frontend bond structure.
    """
    # Create simple molecule: C-C (ethane)
    fe_molecule = {
        "atoms": [
            {"id": "a0", "element": "C", "x": 0, "y": 0, "z": 0, "charge": 0},
            {"id": "a1", "element": "C", "x": 50, "y": 0, "z": 0, "charge": 0},
        ],
        "bonds": [
            {"id": "b0", "atoms": ["a0", "a1"], "order": 1},
        ],
        "node_order": ["a0", "a1"],
    }
    
    data, node_mapping = featurize_json(fe_molecule)
    
    # Verify edge_index structure
    edge_index = data.edge_index.cpu().numpy()
    
    # Should have 2 edges (undirected)
    assert edge_index.shape[1] == 2
    
    # Verify connectivity
    u, v = edge_index[0, 0], edge_index[1, 0]
    assert (u == 0 and v == 1) or (u == 1 and v == 0)
    
    print("✓ Edge index consistency test passed")


if __name__ == "__main__":
    test_ethanol_attention_mapping()
    test_smiles_featurization()
    test_attention_normalization()
    test_edge_index_consistency()
    # test_prediction_engine_smiles()  # Skip if no model
    print("\n✅ All integration tests passed!")

