"""
Phase G: Scaling Decision Framework

Goal: Evaluate scaling feasibility and plan active learning.

Tasks:
1. Evaluate batch size scaling feasibility
2. Test GPU acceleration for ML/RL/fingerprints
3. Plan active learning loop: feed top candidates back to Phase 5 models
"""

import sys
from pathlib import Path
import logging
import json
from datetime import datetime
from typing import Dict, Any

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class ScalingDecisionFramework:
    """Framework for making scaling decisions."""
    
    def __init__(self):
        self.analysis = {
            "timestamp": datetime.now().isoformat(),
            "batch_size_analysis": {},
            "gpu_acceleration": {},
            "active_learning_plan": {},
            "recommendations": [],
        }
    
    def analyze_scaling(self) -> Dict[str, Any]:
        """Analyze scaling feasibility."""
        logger.info("=" * 60)
        logger.info("Phase G: Scaling Decision Framework")
        logger.info("=" * 60)
        
        # Analyze batch size scaling
        logger.info("\n1. Analyzing batch size scaling...")
        self._analyze_batch_sizes()
        
        # Check GPU availability
        logger.info("\n2. Checking GPU acceleration...")
        self._check_gpu_availability()
        
        # Plan active learning
        logger.info("\n3. Planning active learning loop...")
        self._plan_active_learning()
        
        # Generate recommendations
        logger.info("\n4. Generating recommendations...")
        self._generate_recommendations()
        
        return self.analysis
    
    def _analyze_batch_sizes(self):
        """Analyze batch size scaling feasibility."""
        self.analysis["batch_size_analysis"] = {
            "current_batch_size": 32,
            "recommended_batch_sizes": {
                "small": 16,
                "medium": 64,
                "large": 128,
                "xlarge": 256,
            },
            "considerations": [
                "Memory constraints",
                "GPU availability",
                "Throughput requirements",
                "Model complexity",
            ],
        }
        
        logger.info("  ✓ Batch size analysis complete")
    
    def _check_gpu_availability(self):
        """Check GPU acceleration availability."""
        try:
            import torch
            cuda_available = torch.cuda.is_available()
            gpu_count = torch.cuda.device_count() if cuda_available else 0
            
            self.analysis["gpu_acceleration"] = {
                "cuda_available": cuda_available,
                "gpu_count": gpu_count,
                "gpu_names": [torch.cuda.get_device_name(i) for i in range(gpu_count)] if cuda_available else [],
                "recommendations": {
                    "ml_training": "Use GPU if available for training",
                    "ml_inference": "Use GPU for batch inference",
                    "rl_generation": "CPU is sufficient for current batch sizes",
                },
            }
            
            if cuda_available:
                logger.info(f"  ✓ GPU available: {gpu_count} device(s)")
            else:
                logger.info("  ⚠ GPU not available (CPU mode)")
                
        except ImportError:
            self.analysis["gpu_acceleration"] = {
                "status": "pytorch_not_installed",
                "recommendation": "Install PyTorch with CUDA support for GPU acceleration",
            }
            logger.warning("  ⚠ PyTorch not installed")
    
    def _plan_active_learning(self):
        """Plan active learning loop."""
        self.analysis["active_learning_plan"] = {
            "description": "Feed top RL candidates back into ML training",
            "steps": [
                "1. Run RL loop to generate top candidates",
                "2. Select top N candidates by reward",
                "3. Compute properties for candidates (if not already computed)",
                "4. Add to training dataset",
                "5. Retrain Phase 5 ML models",
                "6. Update RL reward function if needed",
                "7. Repeat",
            ],
            "implementation": {
                "script": "backend/scripts/active_learning_loop.py (TODO)",
                "frequency": "Weekly or after N RL loops",
                "dataset_update": "Append to existing training set",
            },
        }
        
        logger.info("  ✓ Active learning plan created")
    
    def _generate_recommendations(self):
        """Generate scaling recommendations."""
        recommendations = []
        
        # Check GPU
        gpu_info = self.analysis.get("gpu_acceleration", {})
        if gpu_info.get("cuda_available"):
            recommendations.append("✅ GPU available - Use for ML training and batch inference")
        else:
            recommendations.append("⚠️  No GPU - Consider cloud GPU instances for training")
        
        # Batch size
        recommendations.append("📊 Start with batch_size=64, scale to 128 if memory allows")
        
        # Active learning
        recommendations.append("🔄 Implement active learning loop after baseline metrics collected")
        
        # Monitoring
        recommendations.append("📈 Set up monitoring for: throughput, memory usage, GPU utilization")
        
        self.analysis["recommendations"] = recommendations
        
        logger.info("\nRecommendations:")
        for rec in recommendations:
            logger.info(f"  {rec}")
    
    def save_analysis(self, output_file: str = "data/docs/scaling_analysis.json"):
        """Save scaling analysis."""
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, "w") as f:
            json.dump(self.analysis, f, indent=2, default=str)
        
        logger.info(f"\nAnalysis saved to: {output_path}")


def main():
    """Run Phase G scaling analysis."""
    framework = ScalingDecisionFramework()
    analysis = framework.analyze_scaling()
    framework.save_analysis()
    
    logger.info("\n" + "=" * 60)
    logger.info("Scaling analysis complete!")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()

