"""
Orchestrator - Main orchestration engine

Coordinates agents and manages workflows.
"""

from typing import Dict, List, Optional, Any
import logging
import time
import uuid

from .agent_registry import AgentRegistry
from .task_router import TaskRouter
from .agent_protocols import Task, TaskResult, Message

logger = logging.getLogger(__name__)


class Orchestrator:
    """
    Main orchestrator for multi-agent drug discovery workflows.
    
    Manages agents, routes tasks, and coordinates workflows.
    """
    
    def __init__(self):
        self.registry = AgentRegistry()
        self.router = TaskRouter(self.registry)
        self.active_tasks: Dict[str, Task] = {}
        self.task_history: List[Dict[str, Any]] = []
    
    def register_agent(self, agent):
        """Register an agent."""
        self.registry.register(agent)
    
    def submit_task(
        self,
        task_type: str,
        input_data: Dict[str, Any],
        config: Optional[Dict[str, Any]] = None,
        priority: int = 0
    ) -> str:
        """
        Submit a task for execution.
        
        Args:
            task_type: Type of task
            input_data: Task input data
            config: Task configuration
            priority: Task priority (higher = more important)
        
        Returns:
            Task ID
        """
        task_id = str(uuid.uuid4())
        
        task = Task(
            task_id=task_id,
            task_type=task_type,
            input_data=input_data,
            config=config or {},
            priority=priority,
        )
        
        self.active_tasks[task_id] = task
        
        logger.info(f"Submitted task: {task_id} (type: {task_type})")
        
        return task_id
    
    async def execute_task(
        self,
        task_id: str,
        routing_strategy: str = "round_robin"
    ) -> TaskResult:
        """
        Execute a task.
        
        Args:
            task_id: Task ID
            routing_strategy: Routing strategy
        
        Returns:
            TaskResult
        """
        task = self.active_tasks.get(task_id)
        if not task:
            return TaskResult(
                task_id=task_id,
                success=False,
                error=f"Task not found: {task_id}",
            )
        
        # Route task to agent
        agent = self.router.route(task, strategy=routing_strategy)
        if not agent:
            result = TaskResult(
                task_id=task_id,
                success=False,
                error=f"No agent available for task type: {task.task_type}",
            )
            self._record_task_result(task, result)
            return result
        
        # Execute task
        start_time = time.time()
        try:
            logger.info(f"Executing task {task_id} with agent {agent.get_name()}")
            result = await agent.run(task)
            execution_time = time.time() - start_time
            result.execution_time = execution_time
        except Exception as e:
            logger.error(f"Error executing task {task_id}: {e}", exc_info=True)
            result = TaskResult(
                task_id=task_id,
                success=False,
                error=str(e),
                execution_time=time.time() - start_time,
            )
        
        # Record result
        self._record_task_result(task, result)
        
        # Remove from active tasks
        if task_id in self.active_tasks:
            del self.active_tasks[task_id]
        
        return result
    
    async def execute_workflow(
        self,
        workflow: List[Dict[str, Any]],
        routing_strategy: str = "round_robin"
    ) -> List[TaskResult]:
        """
        Execute a workflow (sequence of tasks).
        
        Args:
            workflow: List of task definitions
            routing_strategy: Routing strategy
        
        Returns:
            List of TaskResults
        """
        results = []
        
        for task_def in workflow:
            task_type = task_def.get('task_type')
            input_data = task_def.get('input_data', {})
            config = task_def.get('config', {})
            priority = task_def.get('priority', 0)
            
            # Submit task
            task_id = self.submit_task(
                task_type=task_type,
                input_data=input_data,
                config=config,
                priority=priority
            )
            
            # Execute task
            result = await self.execute_task(task_id, routing_strategy=routing_strategy)
            results.append(result)
            
            # Stop workflow if task failed
            if not result.success:
                logger.warning(f"Workflow stopped due to task failure: {task_id}")
                break
        
        return results
    
    def _record_task_result(self, task: Task, result: TaskResult):
        """Record task result in history."""
        self.task_history.append({
            'task_id': task.task_id,
            'task_type': task.task_type,
            'success': result.success,
            'execution_time': result.execution_time,
            'timestamp': time.time(),
        })
        
        # Limit history size
        if len(self.task_history) > 1000:
            self.task_history = self.task_history[-1000:]
    
    def get_task_status(self, task_id: str) -> Optional[Dict[str, Any]]:
        """Get status of a task."""
        if task_id in self.active_tasks:
            return {
                'task_id': task_id,
                'status': 'pending',
                'task_type': self.active_tasks[task_id].task_type,
            }
        
        # Check history
        for record in self.task_history:
            if record['task_id'] == task_id:
                return {
                    'task_id': task_id,
                    'status': 'completed' if record['success'] else 'failed',
                    'task_type': record['task_type'],
                    'execution_time': record['execution_time'],
                    'timestamp': record['timestamp'],
                }
        
        return None
    
    def get_stats(self) -> Dict:
        """Get orchestrator statistics."""
        return {
            'active_tasks': len(self.active_tasks),
            'task_history_size': len(self.task_history),
            'registry_stats': self.registry.get_stats(),
            'routing_stats': self.router.get_routing_stats(),
        }
