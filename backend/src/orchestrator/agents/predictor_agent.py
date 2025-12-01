"""
Predictor Agent - ML prediction agent

Calls ML engine from Phase 5 for property prediction.
"""

import time
import logging
from typing import Dict, Any

from ..agent_protocols import Agent, Task, TaskResult

logger = logging.getLogger(__name__)


class PredictorAgent(Agent):
    """
    Agent for molecular property prediction.
    
    Integrates with ML engine from Phase 5.
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        super().__init__("predictor", config or {})
        self.config.setdefault('capabilities', ['predict', 'property_prediction'])
        self.config.setdefault('model_id', None)
    
    async def run(self, task: Task) -> TaskResult:
        """
        Execute prediction task.
        
        Expected input_data:
            - smiles: SMILES string
            - properties: Optional list of properties to predict
        """
        start_time = time.time()
        
        try:
            input_data = task.input_data
            smiles = input_data.get('smiles')
            
            if not smiles:
                return self._create_result(
                    task.task_id,
                    success=False,
                    error="Missing 'smiles' in input_data",
                )
            
            # Import ML engine from Phase 5
            import sys
            from pathlib import Path
            backend_path = Path(__file__).parent.parent.parent.parent
            if str(backend_path) not in sys.path:
                sys.path.insert(0, str(backend_path))
            
            from ml.prediction_engine import PredictionEngine
            from ml.registry import ModelRegistry
            
            # Initialize engine
            registry = ModelRegistry()
            engine = PredictionEngine(registry)
            
            # Get model ID from config or task
            model_id = task.config.get('model_id') or self.config.get('model_id')
            properties = input_data.get('properties')
            
            # Run prediction
            result = engine.predict(
                input_data={'smiles': smiles},
                model_id=model_id,
                properties=properties,
                return_attention=False,
            )
            
            execution_time = time.time() - start_time
            
            return self._create_result(
                task.task_id,
                success=True,
                output_data={
                    'predictions': result.predictions,
                    'model_id': result.model_id,
                    'confidence': result.confidence,
                },
                execution_time=execution_time,
                metadata=result.metadata or {},
            )
        
        except Exception as e:
            logger.error(f"Prediction agent error: {e}", exc_info=True)
            return self._create_result(
                task.task_id,
                success=False,
                error=str(e),
                execution_time=time.time() - start_time,
            )
