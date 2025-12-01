"""
MD Agent - Molecular dynamics agent

Calls Phase 8 MD engine (mock).
"""

import time
import logging
from typing import Dict, Any

from ..agent_protocols import Agent, Task, TaskResult

logger = logging.getLogger(__name__)


class MDAgent(Agent):
    """
    Agent for molecular dynamics simulations.
    
    Integrates with Phase 8 MD engine.
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        super().__init__("md", config or {})
        self.config.setdefault('capabilities', ['md_simulate'])
        self.config.setdefault('default_steps', 1000)
        self.config.setdefault('default_timestep', 0.001)
        self.config.setdefault('default_temperature', 300.0)
        self.config.setdefault('default_forcefield', 'UFF')
    
    async def run(self, task: Task) -> TaskResult:
        """
        Execute MD simulation task.
        
        Expected input_data:
            - smiles: SMILES string
            - coordinates: Optional initial 3D coordinates
            - steps: Number of MD steps (optional)
            - timestep: Time step in picoseconds (optional)
            - temperature: Temperature in Kelvin (optional)
            - forcefield: Force field name (optional)
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
            
            # Import Phase 8 MD engine
            import sys
            from pathlib import Path
            backend_path = Path(__file__).parent.parent.parent.parent
            if str(backend_path) not in sys.path:
                sys.path.insert(0, str(backend_path))
            
            from src.md import MDEngine
            
            engine = MDEngine()
            
            # Get parameters from input or defaults
            coordinates = input_data.get('coordinates')
            steps = input_data.get('steps', self.config.get('default_steps'))
            timestep = input_data.get('timestep', self.config.get('default_timestep'))
            temperature = input_data.get('temperature', self.config.get('default_temperature'))
            forcefield = input_data.get('forcefield', self.config.get('default_forcefield'))
            
            # Run MD simulation
            result = await engine.simulate(
                smiles=smiles,
                coordinates=coordinates,
                steps=steps,
                timestep=timestep,
                temperature=temperature,
                forcefield=forcefield
            )
            
            # Convert trajectory frames to dicts
            trajectory = [
                {
                    'step': frame.step,
                    'time': frame.time,
                    'coordinates': frame.coordinates,
                    'energy': frame.energy,
                    'temperature': frame.temperature,
                }
                for frame in result.trajectory
            ]
            
            execution_time = time.time() - start_time
            
            return self._create_result(
                task.task_id,
                success=True,
                output_data={
                    'trajectory': trajectory,
                    'final_coordinates': result.final_coordinates,
                    'total_energy': result.total_energy,
                    'steps': steps,
                    'timestep': timestep,
                    'temperature': temperature,
                    'forcefield': forcefield,
                },
                execution_time=execution_time,
                metadata=result.metadata or {},
            )
        
        except Exception as e:
            logger.error(f"MD agent error: {e}", exc_info=True)
            return self._create_result(
                task.task_id,
                success=False,
                error=str(e),
                execution_time=time.time() - start_time,
            )

