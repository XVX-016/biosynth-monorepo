"""
QM Agent - Quantum chemistry agent

Calls Phase 8 QM engine (mock).
"""

import time
import logging
from typing import Dict, Any

from ..agent_protocols import Agent, Task, TaskResult

logger = logging.getLogger(__name__)


class QMAgent(Agent):
    """
    Agent for quantum chemistry calculations.
    
    Integrates with Phase 8 QM engine.
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        super().__init__("qm", config or {})
        self.config.setdefault('capabilities', ['qm_energy', 'qm_optimize', 'qm_properties'])
        self.config.setdefault('default_method', 'B3LYP')
        self.config.setdefault('default_basis', '6-31G')
    
    async def run(self, task: Task) -> TaskResult:
        """
        Execute QM task.
        
        Expected input_data:
            - smiles: SMILES string
            - task_subtype: 'energy', 'optimize', or 'properties'
            - coordinates: Optional 3D coordinates
            - method: QM method (optional)
            - basis: Basis set (optional)
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
            
            # Import Phase 8 QM engine
            import sys
            from pathlib import Path
            backend_path = Path(__file__).parent.parent.parent.parent
            if str(backend_path) not in sys.path:
                sys.path.insert(0, str(backend_path))
            
            from src.qm import QMEngine
            
            engine = QMEngine()
            
            task_subtype = input_data.get('task_subtype', 'energy')
            method = input_data.get('method', self.config.get('default_method'))
            basis = input_data.get('basis', self.config.get('default_basis'))
            coordinates = input_data.get('coordinates')
            
            if task_subtype == 'energy':
                result = await engine.compute_energy(
                    smiles=smiles,
                    coordinates=coordinates,
                    method=method,
                    basis=basis
                )
                output_data = {
                    'energy': result.energy,
                    'coordinates': result.coordinates,
                    'method': method,
                    'basis': basis,
                }
            
            elif task_subtype == 'optimize':
                result = await engine.optimize_geometry(
                    smiles=smiles,
                    initial_coordinates=coordinates,
                    method=method,
                    basis=basis
                )
                output_data = {
                    'energy': result.energy,
                    'coordinates': result.coordinates,
                    'forces': result.forces,
                    'method': method,
                    'basis': basis,
                }
            
            elif task_subtype == 'properties':
                result = await engine.compute_properties(
                    smiles=smiles,
                    coordinates=coordinates,
                    method=method,
                    basis=basis
                )
                output_data = {
                    'energy': result.energy,
                    'dipole': result.dipole,
                    'charges': result.charges,
                    'method': method,
                    'basis': basis,
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
                metadata=result.metadata or {},
            )
        
        except Exception as e:
            logger.error(f"QM agent error: {e}", exc_info=True)
            return self._create_result(
                task.task_id,
                success=False,
                error=str(e),
                execution_time=time.time() - start_time,
            )

