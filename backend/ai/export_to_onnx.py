"""
Export PropertyPredictor model to ONNX format
"""
import os
import torch
import onnx
import onnxruntime as ort
from backend.models.property_predictor import create_model


def export_to_onnx(
    weights_path: str = 'backend/weights/property_predictor.pt',
    output_path: str = 'backend/weights/predictor.onnx',
    input_dim: int = 2048
):
    """
    Export trained PyTorch model to ONNX format
    
    Args:
        weights_path: Path to PyTorch model weights
        output_path: Path to save ONNX model
        input_dim: Dimension of input fingerprint
    """
    # Load PyTorch model
    model = create_model(input_dim=input_dim)
    
    if os.path.exists(weights_path):
        model.load_state_dict(torch.load(weights_path, map_location='cpu'))
        print(f"Loaded weights from {weights_path}")
    else:
        print(f"Warning: Weights not found at {weights_path}. Exporting untrained model.")
    
    model.eval()
    
    # Create dummy input
    dummy_input = torch.randn(1, input_dim)
    
    # Export to ONNX
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    torch.onnx.export(
        model,
        dummy_input,
        output_path,
        input_names=['fingerprint'],
        output_names=['properties'],
        dynamic_axes={
            'fingerprint': {0: 'batch_size'},
            'properties': {0: 'batch_size'}
        },
        opset_version=11,
        do_constant_folding=True
    )
    
    print(f"Model exported to {output_path}")
    
    # Verify ONNX model
    onnx_model = onnx.load(output_path)
    onnx.checker.check_model(onnx_model)
    print("ONNX model verification passed")
    
    # Test inference
    session = ort.InferenceSession(output_path)
    input_name = session.get_inputs()[0].name
    output_name = session.get_outputs()[0].name
    
    # Run test inference
    test_input = dummy_input.numpy()
    result = session.run([output_name], {input_name: test_input})
    print(f"Test inference successful. Output shape: {result[0].shape}")
    
    return output_path


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Export PropertyPredictor to ONNX')
    parser.add_argument('--weights', type=str, default='backend/weights/property_predictor.pt',
                       help='Path to PyTorch model weights')
    parser.add_argument('--output', type=str, default='backend/weights/predictor.onnx',
                       help='Path to save ONNX model')
    parser.add_argument('--input-dim', type=int, default=2048,
                       help='Input dimension (fingerprint size)')
    
    args = parser.parse_args()
    
    export_to_onnx(
        weights_path=args.weights,
        output_path=args.output,
        input_dim=args.input_dim
    )

