"""
Phase F: Documentation & Reproducibility

Goal: Prevent lost knowledge and debugging chaos.

Tasks:
1. Document every pipeline step, input/output, and API
2. Record environment, dependencies, and configurations
3. Save sample datasets and test runs
"""

import sys
from pathlib import Path
import logging
import json
import subprocess
from datetime import datetime
from typing import Dict, Any

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class DocumentationGenerator:
    """Generates comprehensive documentation for reproducibility."""
    
    def __init__(self):
        self.docs = {
            "timestamp": datetime.now().isoformat(),
            "environment": {},
            "dependencies": {},
            "pipeline_steps": {},
            "api_endpoints": {},
            "sample_datasets": {},
        }
    
    def generate_documentation(self) -> Dict[str, Any]:
        """Generate all documentation."""
        logger.info("=" * 60)
        logger.info("Phase F: Documentation & Reproducibility")
        logger.info("=" * 60)
        
        # Collect environment info
        logger.info("\n1. Collecting environment information...")
        self._collect_environment()
        
        # Collect dependencies
        logger.info("\n2. Collecting dependencies...")
        self._collect_dependencies()
        
        # Document pipeline steps
        logger.info("\n3. Documenting pipeline steps...")
        self._document_pipeline_steps()
        
        # Document API endpoints
        logger.info("\n4. Documenting API endpoints...")
        self._document_api_endpoints()
        
        # Document sample datasets
        logger.info("\n5. Documenting sample datasets...")
        self._document_sample_datasets()
        
        return self.docs
    
    def _collect_environment(self):
        """Collect environment information."""
        import platform
        import sys
        import os
        
        self.docs["environment"] = {
            "platform": platform.platform(),
            "python_version": sys.version,
            "python_executable": sys.executable,
            "working_directory": str(Path.cwd()),
            "environment_variables": {
                "PATH": os.environ.get("PATH", "")[:100] + "...",  # Truncate
            },
        }
        
        logger.info("  ✓ Environment info collected")
    
    def _collect_dependencies(self):
        """Collect dependency information."""
        try:
            # Try to get pip list
            result = subprocess.run(
                ["pip", "list", "--format=json"],
                capture_output=True,
                text=True,
                timeout=10,
            )
            if result.returncode == 0:
                dependencies = json.loads(result.stdout)
                self.docs["dependencies"] = {
                    "python_packages": dependencies[:20],  # First 20
                    "total_count": len(dependencies),
                }
                logger.info(f"  ✓ Found {len(dependencies)} Python packages")
            else:
                self.docs["dependencies"] = {"status": "failed", "error": result.stderr}
        except Exception as e:
            self.docs["dependencies"] = {"status": "error", "error": str(e)}
            logger.warning(f"  ⚠ Could not collect dependencies: {e}")
    
    def _document_pipeline_steps(self):
        """Document pipeline steps."""
        self.docs["pipeline_steps"] = {
            "dataset_preparation": {
                "script": "backend/scripts/prepare_datasets.py",
                "steps": [
                    "collect_property_datasets",
                    "prepare_screening_libraries",
                    "standardize_molecules",
                    "split_datasets",
                    "generate_fingerprints",
                    "generate_mock_rewards",
                ],
            },
            "ml_training": {
                "script": "backend/scripts/train_ml_models.py",
                "steps": [
                    "load_datasets",
                    "configure_models",
                    "train_models",
                    "evaluate_models",
                    "register_models",
                ],
            },
            "rl_loop": {
                "script": "backend/scripts/run_rl_loop_test.py",
                "description": "Run small RL loop test",
            },
        }
        
        logger.info("  ✓ Pipeline steps documented")
    
    def _document_api_endpoints(self):
        """Document API endpoints."""
        self.docs["api_endpoints"] = {
            "phase10": {
                "generate": "POST /api/phase10/generate",
                "evaluate": "POST /api/phase10/evaluate",
                "run_loop": "POST /api/phase10/run_loop",
                "top_candidates": "GET /api/phase10/top_candidates",
                "statistics": "GET /api/phase10/statistics",
                "iteration_logs": "GET /api/phase10/iteration_logs",
            },
            "ml": {
                "predict": "POST /api/predict/property",
                "attention_map": "POST /api/predict/attention-map",
            },
            "search": {
                "similarity": "GET /api/search/similarity",
                "substructure": "GET /api/search/substructure",
            },
            "conformers": {
                "generate": "POST /api/conformers/generate",
            },
        }
        
        logger.info("  ✓ API endpoints documented")
    
    def _document_sample_datasets(self):
        """Document sample datasets."""
        data_dir = Path("data")
        
        datasets = {}
        if data_dir.exists():
            for subdir in ["datasets", "libraries", "models"]:
                subdir_path = data_dir / subdir
                if subdir_path.exists():
                    files = list(subdir_path.rglob("*"))
                    datasets[subdir] = {
                        "path": str(subdir_path),
                        "file_count": len([f for f in files if f.is_file()]),
                    }
        
        self.docs["sample_datasets"] = datasets
        
        logger.info("  ✓ Sample datasets documented")
    
    def save_documentation(self, output_file: str = "data/docs/reproducibility_docs.json"):
        """Save documentation."""
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, "w") as f:
            json.dump(self.docs, f, indent=2, default=str)
        
        logger.info(f"\nDocumentation saved to: {output_path}")


def main():
    """Run Phase F documentation generation."""
    generator = DocumentationGenerator()
    docs = generator.generate_documentation()
    generator.save_documentation()
    
    logger.info("\n" + "=" * 60)
    logger.info("Documentation generation complete!")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()

