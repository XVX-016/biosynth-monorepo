"""
Phase 9: Multi-Agent Orchestrator API routes

Provides endpoints for task submission, workflow execution, and agent management.
"""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field
from typing import List, Optional, Dict, Any
import logging
import sys
from pathlib import Path

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

logger = logging.getLogger(__name__)

router = APIRouter()

# Global orchestrator instance
_orchestrator = Orchestrator()

# Register default agents
_orchestrator.register_agent(PredictorAgent())
_orchestrator.register_agent(ScreeningAgent())
_orchestrator.register_agent(QMAgent())
_orchestrator.register_agent(MDAgent())


# Request models
class TaskSubmitRequest(BaseModel):
    """Request to submit a task"""
    task_type: str = Field(..., description="Type of task")
    input_data: Dict[str, Any] = Field(..., description="Task input data")
    config: Dict[str, Any] = Field(default_factory=dict, description="Task configuration")
    priority: int = Field(0, description="Task priority")


class WorkflowExecuteRequest(BaseModel):
    """Request to execute a workflow"""
    workflow: List[Dict[str, Any]] = Field(..., description="Workflow definition")
    routing_strategy: str = Field("round_robin", description="Routing strategy")


@router.post("/task/submit")
async def submit_task(request: TaskSubmitRequest):
    """
    Submit a task for execution.
    """
    try:
        task_id = _orchestrator.submit_task(
            task_type=request.task_type,
            input_data=request.input_data,
            config=request.config,
            priority=request.priority
        )
        
        return {
            "task_id": task_id,
            "status": "submitted",
        }
    except Exception as e:
        logger.error(f"Task submission error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/task/execute")
async def execute_task(
    task_id: str,
    routing_strategy: str = "round_robin"
):
    """
    Execute a submitted task.
    """
    try:
        result = await _orchestrator.execute_task(task_id, routing_strategy=routing_strategy)
        
        return {
            "task_id": result.task_id,
            "success": result.success,
            "output_data": result.output_data,
            "error": result.error,
            "execution_time": result.execution_time,
            "metadata": result.metadata,
        }
    except Exception as e:
        logger.error(f"Task execution error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/task/status/{task_id}")
async def get_task_status(task_id: str):
    """
    Get status of a task.
    """
    try:
        status = _orchestrator.get_task_status(task_id)
        if not status:
            raise HTTPException(status_code=404, detail=f"Task not found: {task_id}")
        return status
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Get task status error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/workflow/execute")
async def execute_workflow(request: WorkflowExecuteRequest):
    """
    Execute a workflow (sequence of tasks).
    """
    try:
        results = await _orchestrator.execute_workflow(
            workflow=request.workflow,
            routing_strategy=request.routing_strategy
        )
        
        return {
            "results": [
                {
                    "task_id": r.task_id,
                    "success": r.success,
                    "output_data": r.output_data,
                    "error": r.error,
                    "execution_time": r.execution_time,
                }
                for r in results
            ],
            "count": len(results),
            "all_successful": all(r.success for r in results),
        }
    except Exception as e:
        logger.error(f"Workflow execution error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/agents")
async def list_agents():
    """
    List all registered agents.
    """
    try:
        stats = _orchestrator.registry.get_stats()
        return {
            "agents": stats['agents'],
            "num_agents": stats['num_agents'],
            "task_routing": stats['task_routing'],
        }
    except Exception as e:
        logger.error(f"List agents error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/stats")
async def get_orchestrator_stats():
    """
    Get orchestrator statistics.
    """
    try:
        return _orchestrator.get_stats()
    except Exception as e:
        logger.error(f"Get stats error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))

