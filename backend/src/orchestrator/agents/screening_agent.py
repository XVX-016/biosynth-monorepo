"""
Screening Agent - Molecular screening agent

Calls Phase 7 screening engine.
"""

import time
import logging
from typing import Dict, Any

from ..agent_protocols import Agent, Task, TaskResult

logger = logging.getLogger(__name__)


class ScreeningAgent(Agent):
    """
    Agent for molecular screening.
    
    Integrates with Phase 7 screening engine.
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        super().__init__("screening", config or {})
        self.config.setdefault('capabilities', ['screen', 'similarity_search', 'substructure_search'])
    
    async def run(self, task: Task) -> TaskResult:
        """
        Execute screening task.
        
        Expected input_data:
            - query_smiles: Query SMILES string
            - task_subtype: 'similarity' or 'substructure'
            - k: Number of results (for similarity)
            - threshold: Similarity threshold (for similarity)
            - smarts: SMARTS pattern (for substructure)
        """
        start_time = time.time()
        
        try:
            input_data = task.input_data
            task_subtype = input_data.get('task_subtype', 'similarity')
            
            # Import Phase 7 search engine
            import sys
            from pathlib import Path
            backend_path = Path(__file__).parent.parent.parent.parent
            if str(backend_path) not in sys.path:
                sys.path.insert(0, str(backend_path))
            
            from src.search import SearchEngine
            
            engine = SearchEngine()
            
            if task_subtype == 'similarity':
                query_smiles = input_data.get('query_smiles')
                if not query_smiles:
                    return self._create_result(
                        task.task_id,
                        success=False,
                        error="Missing 'query_smiles' in input_data",
                    )
                
                k = input_data.get('k', 10)
                threshold = input_data.get('threshold')
                
                results = engine.similarity_search(
                    smiles=query_smiles,
                    k=k,
                    threshold=threshold
                )
                
                output_data = {
                    'results': results,
                    'count': len(results),
                    'query_smiles': query_smiles,
                }
            
            elif task_subtype == 'substructure':
                smarts = input_data.get('smarts')
                if not smarts:
                    return self._create_result(
                        task.task_id,
                        success=False,
                        error="Missing 'smarts' in input_data",
                    )
                
                max_results = input_data.get('max_results', 100)
                
                results = engine.substructure_search(
                    smarts=smarts,
                    max_results=max_results
                )
                
                output_data = {
                    'results': results,
                    'count': len(results),
                    'query_smarts': smarts,
                }
            
            else:
                return self._create_result(
                    task.task_id,
                    success=False,
                    error=f"Unknown task_subtype: {task_subtype}",
                )
            
            execution_time = time.time() - start_time
            
            return self._create_result(
                task.task_id,
                success=True,
                output_data=output_data,
                execution_time=execution_time,
            )
        
        except Exception as e:
            logger.error(f"Screening agent error: {e}", exc_info=True)
            return self._create_result(
                task.task_id,
                success=False,
                error=str(e),
                execution_time=time.time() - start_time,
            )
