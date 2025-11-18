"""
Test ONNX predictor
"""
import pytest
import numpy as np
import os
from backend.ai.onnx_predictor import ONNXPredictor


def test_onnx_predictor_requires_file():
    """Test that ONNX predictor requires model file"""
    with pytest.raises(FileNotFoundError):
        ONNXPredictor(onnx_path='nonexistent.onnx')


@pytest.mark.skipif(
    not os.path.exists('backend/weights/predictor.onnx'),
    reason="ONNX model not available"
)
def test_onnx_predictor_loads():
    """Test that ONNX predictor loads model"""
    predictor = ONNXPredictor()
    assert predictor.session is not None
    assert predictor.input_name is not None
    assert predictor.output_name is not None


@pytest.mark.skipif(
    not os.path.exists('backend/weights/predictor.onnx'),
    reason="ONNX model not available"
)
def test_onnx_predictor_predict():
    """Test ONNX predictor makes predictions"""
    predictor = ONNXPredictor()
    
    # Create dummy features
    features = [0.0] * 2048
    features[0] = 1.0
    
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


@pytest.mark.skipif(
    not os.path.exists('backend/weights/predictor.onnx'),
    reason="ONNX model not available"
)
def test_onnx_matches_pytorch():
    """Test that ONNX predictions match PyTorch (within tolerance)"""
    from backend.ai.predictor import PropertyPredictor
    
    # Create features
    features = [0.5] * 2048
    
    # PyTorch prediction
    pytorch_predictor = PropertyPredictor()
    pytorch_result = pytorch_predictor.predict(features)
    
    # ONNX prediction
    onnx_predictor = ONNXPredictor()
    onnx_result = onnx_predictor.predict(features)
    
    # Compare (allow small tolerance for numerical differences)
    tolerance = 0.01
    assert abs(pytorch_result['stability'] - onnx_result['stability']) < tolerance
    assert abs(pytorch_result['toxicity'] - onnx_result['toxicity']) < tolerance
    assert abs(pytorch_result['solubility'] - onnx_result['solubility']) < tolerance

