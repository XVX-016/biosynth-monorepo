"""
Task Router - Routes tasks to appropriate agents

Supports round-robin and rule-based routing strategies.
"""

from typing import List, Optional, Dict, Any
import logging
import random

from .agent_protocols import Task
from .agent_registry import AgentRegistry
from .agent_protocols import Agent

logger = logging.getLogger(__name__)


class TaskRouter:
    """
    Routes tasks to agents based on strategy.
    
    Supports:
    - Round-robin: Distribute tasks evenly
    - Rule-based: Route based on task properties
    - Load-based: Route to least loaded agent (future)
    """
    
    def __init__(self, registry: AgentRegistry):
        """
        Args:
            registry: AgentRegistry instance
        """
        self.registry = registry
        self.round_robin_counters: Dict[str, int] = {}  # task_type -> counter
    
    def route(
        self,
        task: Task,
        strategy: str = "round_robin"
    ) -> Optional[Agent]:
        """
        Route task to an agent.
        
        Args:
            task: Task to route
            strategy: Routing strategy ('round_robin', 'rule_based', 'random')
        
        Returns:
            Agent instance or None if no agent available
        """
        # Get agents that can handle this task type
        agents = self.registry.get_agents_for_task(task.task_type)
        
        if not agents:
            logger.warning(f"No agents available for task type: {task.task_type}")
            return None
        
        if len(agents) == 1:
            return agents[0]
        
        # Apply routing strategy
        if strategy == "round_robin":
            return self._round_robin_route(task.task_type, agents)
        elif strategy == "rule_based":
            return self._rule_based_route(task, agents)
        elif strategy == "random":
            return random.choice(agents)
        else:
            logger.warning(f"Unknown routing strategy: {strategy}, using round_robin")
            return self._round_robin_route(task.task_type, agents)
    
    def _round_robin_route(
        self,
        task_type: str,
        agents: List[Agent]
    ) -> Agent:
        """
        Round-robin routing: distribute tasks evenly.
        
        Args:
            task_type: Type of task
            agents: List of available agents
        
        Returns:
            Selected agent
        """
        if task_type not in self.round_robin_counters:
            self.round_robin_counters[task_type] = 0
        
        counter = self.round_robin_counters[task_type]
        selected = agents[counter % len(agents)]
        
        self.round_robin_counters[task_type] = counter + 1
        
        return selected
    
    def _rule_based_route(
        self,
        task: Task,
        agents: List[Agent]
    ) -> Agent:
        """
        Rule-based routing: route based on task properties.
        
        Args:
            task: Task to route
            agents: List of available agents
        
        Returns:
            Selected agent
        """
        # Check task config for preferred agent
        preferred_agent_name = task.config.get('preferred_agent')
        if preferred_agent_name:
            for agent in agents:
                if agent.get_name() == preferred_agent_name:
                    return agent
        
        # Check task priority
        if task.priority > 5:
            # High priority: use first available agent
            return agents[0]
        
        # Default: round-robin
        return self._round_robin_route(task.task_type, agents)
    
    def get_routing_stats(self) -> Dict:
        """Get routing statistics."""
        return {
            'round_robin_counters': self.round_robin_counters.copy(),
        }
