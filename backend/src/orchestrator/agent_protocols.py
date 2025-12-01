"""
Agent Protocols - Base classes and message protocol for agents

Defines the interface that all agents must implement.
"""

from typing import Protocol, Dict, Any, Optional, List
from dataclasses import dataclass, field
from abc import ABC, abstractmethod
import json
import logging

logger = logging.getLogger(__name__)


@dataclass
class Message:
    """
    JSON message protocol for agent communication.
    
    All messages between agents and orchestrator use this format.
    """
    message_id: str
    message_type: str  # 'task', 'result', 'error', 'status'
    sender: str
    receiver: str
    payload: Dict[str, Any] = field(default_factory=dict)
    timestamp: float = field(default_factory=lambda: __import__('time').time())
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def to_json(self) -> str:
        """Serialize message to JSON string."""
        return json.dumps({
            'message_id': self.message_id,
            'message_type': self.message_type,
            'sender': self.sender,
            'receiver': self.receiver,
            'payload': self.payload,
            'timestamp': self.timestamp,
            'metadata': self.metadata,
        })
    
    @classmethod
    def from_json(cls, json_str: str) -> 'Message':
        """Deserialize message from JSON string."""
        data = json.loads(json_str)
        return cls(**data)


@dataclass
class Task:
    """
    Task definition for agents.
    """
    task_id: str
    task_type: str  # 'predict', 'screen', 'qm_optimize', 'md_simulate', etc.
    input_data: Dict[str, Any]
    config: Dict[str, Any] = field(default_factory=dict)
    priority: int = 0  # Higher = more important
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class TaskResult:
    """
    Result from agent task execution.
    """
    task_id: str
    success: bool
    output_data: Dict[str, Any] = field(default_factory=dict)
    error: Optional[str] = None
    execution_time: float = 0.0
    metadata: Dict[str, Any] = field(default_factory=dict)


class AgentProtocol(Protocol):
    """
    Protocol that all agents must implement.
    """
    
    def get_name(self) -> str:
        """Get agent name."""
        ...
    
    def get_config(self) -> Dict[str, Any]:
        """Get agent configuration."""
        ...
    
    async def run(self, task: Task) -> TaskResult:
        """
        Execute a task.
        
        Args:
            task: Task to execute
        
        Returns:
            TaskResult with output data or error
        """
        ...


class Agent(ABC):
    """
    Base class for all agents.
    
    Provides common functionality and enforces interface.
    """
    
    def __init__(self, name: str, config: Optional[Dict[str, Any]] = None):
        """
        Args:
            name: Agent name/identifier
            config: Agent configuration dict
        """
        self.name = name
        self.config = config or {}
        self.logger = logging.getLogger(f"agent.{name}")
    
    def get_name(self) -> str:
        """Get agent name."""
        return self.name
    
    def get_config(self) -> Dict[str, Any]:
        """Get agent configuration."""
        return self.config.copy()
    
    @abstractmethod
    async def run(self, task: Task) -> TaskResult:
        """
        Execute a task.
        
        Must be implemented by subclasses.
        
        Args:
            task: Task to execute
        
        Returns:
            TaskResult with output data or error
        """
        pass
    
    def _create_result(
        self,
        task_id: str,
        success: bool,
        output_data: Optional[Dict[str, Any]] = None,
        error: Optional[str] = None,
        execution_time: float = 0.0,
        metadata: Optional[Dict[str, Any]] = None
    ) -> TaskResult:
        """
        Helper to create TaskResult.
        """
        return TaskResult(
            task_id=task_id,
            success=success,
            output_data=output_data or {},
            error=error,
            execution_time=execution_time,
            metadata=metadata or {},
        )
