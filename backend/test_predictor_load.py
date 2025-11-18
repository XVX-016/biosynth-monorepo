"""
Test model loading
"""
import pytest
import torch
from backend.ai.predictor import ModelLoader, PropertyPredictor
from backend.ai.property_predictor import create_model


def test_model_loader_creates_model():
    """Test that ModelLoader creates a model even without weights"""
    loader = ModelLoader(weights_path='nonexistent.pt')
    assert loader.model is not None
    assert loader.device is not None


def test_property_predictor_initializes():
    """Test that PropertyPredictor initializes"""
    predictor = PropertyPredictor()
    assert predictor.loader is not None


def test_model_forward_pass():
    """Test that model can run forward pass"""
    model = create_model()
    model.eval()
    
    # Create dummy input
    dummy_input = torch.randn(1, 2048)
    
    # Run forward pass
    with torch.no_grad():
        output = model(dummy_input)
    
    # Check output shape
    assert output.shape == (1, 5)  # batch_size=1, output_dim=5


def test_predictor_predict():
    """Test that predictor can make predictions"""
    predictor = PropertyPredictor()
    
    # Create dummy features
    features = [0.0] * 2048
    features[0] = 1.0  # Set one feature
    
    # Make prediction
    result = predictor.predict(features)
    
    # Check result structure
    assert 'stability' in result
    assert 'toxicity' in result
    assert 'solubility' in result
    assert 'bioavailability' in result
    assert 'novelty' in result
    
    # Check types
    assert isinstance(result['stability'], float)
    assert isinstance(result['toxicity'], float)
    assert isinstance(result['solubility'], float)
    assert isinstance(result['bioavailability'], float)
    assert isinstance(result['novelty'], float)

