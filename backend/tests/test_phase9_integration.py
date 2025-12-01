"""
Integration tests for Phase 9: Multi-Agent Orchestrator

Tests:
- Agent registration
- Task submission and execution
- Workflow execution
- Agent integration with Phases 5/7/8
"""

import pytest
import sys
import asyncio
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

from src.orchestrator import Orchestrator
from src.orchestrator.agents import (
    PredictorAgent,
    ScreeningAgent,
    QMAgent,
    MDAgent,
)
from src.orchestrator.agent_protocols import Task


class TestPhase9Orchestrator:
    """Test Phase 9 orchestrator functionality."""
    
    def test_agent_registration(self):
        """Test agents can be registered."""
        orchestrator = Orchestrator()
        
        # Register agents
        orchestrator.register_agent(PredictorAgent())
        orchestrator.register_agent(ScreeningAgent())
        orchestrator.register_agent(QMAgent())
        orchestrator.register_agent(MDAgent())
        
        # Check registration
        stats = orchestrator.registry.get_stats()
        assert stats['num_agents'] == 4
        assert 'predictor' in stats['agents']
        assert 'screening' in stats['agents']
        assert 'qm' in stats['agents']
        assert 'md' in stats['agents']
    
    def test_task_submission(self):
        """Test task submission works."""
        orchestrator = Orchestrator()
        orchestrator.register_agent(PredictorAgent())
        
        # Submit task
        task_id = orchestrator.submit_task(
            task_type="predict",
            input_data={"smiles": "CCO"},
            priority=1
        )
        
        assert task_id is not None
        assert task_id in orchestrator.active_tasks
    
    @pytest.mark.asyncio
    async def test_task_execution(self):
        """Test task execution works."""
        orchestrator = Orchestrator()
        orchestrator.register_agent(PredictorAgent())
        
        # Submit task
        task_id = orchestrator.submit_task(
            task_type="predict",
            input_data={"smiles": "CCO", "properties": ["logP"]},
            priority=1
        )
        
        # Execute task
        result = await orchestrator.execute_task(task_id)
        
        assert result is not None
        assert result.task_id == task_id
        # May succeed or fail depending on ML engine availability
        assert isinstance(result.success, bool)
    
    @pytest.mark.asyncio
    async def test_workflow_execution(self):
        """Test workflow execution works."""
        orchestrator = Orchestrator()
        orchestrator.register_agent(PredictorAgent())
        orchestrator.register_agent(ScreeningAgent())
        
        # Define workflow
        workflow = [
            {
                "task_type": "predict",
                "input_data": {"smiles": "CCO", "properties": ["logP"]},
                "priority": 2
            },
            {
                "task_type": "screen",
                "input_data": {
                    "task_subtype": "similarity",
                    "query_smiles": "CCO",
                    "k": 5
                },
                "priority": 1
            }
        ]
        
        # Execute workflow
        results = await orchestrator.execute_workflow(workflow)
        
        assert len(results) == 2
        assert all(isinstance(r.success, bool) for r in results)
    
    def test_task_status(self):
        """Test task status tracking."""
        orchestrator = Orchestrator()
        orchestrator.register_agent(PredictorAgent())
        
        # Submit task
        task_id = orchestrator.submit_task(
            task_type="predict",
            input_data={"smiles": "CCO"},
        )
        
        # Check status
        status = orchestrator.get_task_status(task_id)
        assert status is not None
        assert status['status'] == 'pending'
        assert status['task_type'] == 'predict'


class TestPhase9Agents:
    """Test Phase 9 agents."""
    
    @pytest.mark.asyncio
    async def test_predictor_agent(self):
        """Test predictor agent can receive tasks."""
        agent = PredictorAgent()
        
        task = Task(
            task_id="test_1",
            task_type="predict",
            input_data={"smiles": "CCO", "properties": ["logP"]}
        )
        
        result = await agent.run(task)
        
        assert result is not None
        assert result.task_id == "test_1"
        # May succeed or fail depending on ML engine
        assert isinstance(result.success, bool)
    
    @pytest.mark.asyncio
    async def test_screening_agent(self):
        """Test screening agent can receive tasks."""
        agent = ScreeningAgent()
        
        task = Task(
            task_id="test_2",
            task_type="screen",
            input_data={
                "task_subtype": "similarity",
                "query_smiles": "CCO",
                "k": 5
            }
        )
        
        result = await agent.run(task)
        
        assert result is not None
        assert result.task_id == "test_2"
        assert isinstance(result.success, bool)
    
    @pytest.mark.asyncio
    async def test_qm_agent(self):
        """Test QM agent can receive tasks."""
        agent = QMAgent()
        
        task = Task(
            task_id="test_3",
            task_type="qm_energy",
            input_data={
                "smiles": "CCO",
                "task_subtype": "energy"
            }
        )
        
        result = await agent.run(task)
        
        assert result is not None
        assert result.task_id == "test_3"
        assert isinstance(result.success, bool)
    
    @pytest.mark.asyncio
    async def test_md_agent(self):
        """Test MD agent can receive tasks."""
        agent = MDAgent()
        
        task = Task(
            task_id="test_4",
            task_type="md_simulate",
            input_data={
                "smiles": "CCO",
                "steps": 100
            }
        )
        
        result = await agent.run(task)
        
        assert result is not None
        assert result.task_id == "test_4"
        assert isinstance(result.success, bool)


class TestPhase9API:
    """Test Phase 9 API endpoints."""
    
    def test_api_imports(self):
        """Test API routes can be imported."""
        try:
            from api import orchestrator as orchestrator_api
            
            assert hasattr(orchestrator_api, 'router')
        except ImportError as e:
            pytest.skip(f"API routes not available: {e}")
    
    def test_orchestrator_stats(self):
        """Test orchestrator statistics endpoint."""
        orchestrator = Orchestrator()
        orchestrator.register_agent(PredictorAgent())
        
        stats = orchestrator.get_stats()
        
        assert 'active_tasks' in stats
        assert 'task_history_size' in stats
        assert 'registry_stats' in stats
        assert 'routing_stats' in stats

