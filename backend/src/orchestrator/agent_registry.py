"""
Agent Registry - Manages registered agents

Tracks available agents and their capabilities.
"""

from typing import Dict, List, Optional
import logging

from .agent_protocols import Agent, AgentProtocol

logger = logging.getLogger(__name__)


class AgentRegistry:
    """
    Registry for managing agents.
    
    Tracks available agents and routes tasks to appropriate agents.
    """
    
    def __init__(self):
        # agent_name -> Agent instance
        self.agents: Dict[str, Agent] = {}
        
        # task_type -> list of agent names that can handle it
        self.task_routing: Dict[str, List[str]] = {}
    
    def register(self, agent: Agent):
        """
        Register an agent.
        
        Args:
            agent: Agent instance to register
        """
        name = agent.get_name()
        self.agents[name] = agent
        
        # Update task routing based on agent capabilities
        config = agent.get_config()
        capabilities = config.get('capabilities', [])
        
        for task_type in capabilities:
            if task_type not in self.task_routing:
                self.task_routing[task_type] = []
            if name not in self.task_routing[task_type]:
                self.task_routing[task_type].append(name)
        
        logger.info(f"Registered agent: {name} with capabilities: {capabilities}")
    
    def unregister(self, agent_name: str):
        """Unregister an agent."""
        if agent_name in self.agents:
            del self.agents[agent_name]
            
            # Remove from task routing
            for task_type, agent_list in self.task_routing.items():
                if agent_name in agent_list:
                    agent_list.remove(agent_name)
            
            logger.info(f"Unregistered agent: {agent_name}")
    
    def get_agent(self, agent_name: str) -> Optional[Agent]:
        """Get agent by name."""
        return self.agents.get(agent_name)
    
    def list_agents(self) -> List[str]:
        """List all registered agent names."""
        return list(self.agents.keys())
    
    def get_agents_for_task(self, task_type: str) -> List[Agent]:
        """
        Get agents that can handle a task type.
        
        Args:
            task_type: Type of task
        
        Returns:
            List of agents that can handle this task
        """
        agent_names = self.task_routing.get(task_type, [])
        return [self.agents[name] for name in agent_names if name in self.agents]
    
    def get_stats(self) -> Dict:
        """Get registry statistics."""
        return {
            'num_agents': len(self.agents),
            'agents': list(self.agents.keys()),
            'task_routing': {
                task_type: len(agents)
                for task_type, agents in self.task_routing.items()
            },
        }
