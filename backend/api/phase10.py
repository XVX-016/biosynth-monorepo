"""
Phase 10 API Endpoints

Endpoints for RL + Generative molecule design:
- /api/phase10/generate
- /api/phase10/evaluate
- /api/phase10/run_loop
- /api/phase10/top_candidates
"""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List, Dict, Optional, Any
import logging
import asyncio

# Import Phase 10 components
from src.phase10 import (
    RLAgent,
    GenerativeAgent,
    RewardFunction,
    Evaluator,
    WorkflowLoop,
    DatasetUtils,
    Phase10Orchestrator,
)

# Import Phase 9 orchestrator for integration
try:
    from src.orchestrator import Orchestrator
    ORCHESTRATOR_AVAILABLE = True
except ImportError:
    ORCHESTRATOR_AVAILABLE = False
    Orchestrator = None

logger = logging.getLogger(__name__)

router = APIRouter()

# Global instances (initialized on startup)
_phase10_orchestrator: Optional[Phase10Orchestrator] = None
_orchestrator: Optional[Any] = None


def init_phase10(orchestrator: Optional[Any] = None):
    """
    Initialize Phase 10 components.
    
    Args:
        orchestrator: Phase 9 Orchestrator instance (optional)
    """
    global _phase10_orchestrator, _orchestrator
    
    _orchestrator = orchestrator
    
    # Initialize components
    rl_agent = RLAgent()
    generative_agent = GenerativeAgent()
    reward_function = RewardFunction()
    dataset_utils = DatasetUtils()
    
    # Initialize evaluator with orchestrator
    evaluator = Evaluator(
        reward_function=reward_function,
        orchestrator=_orchestrator,
    )
    
    # Initialize workflow loop
    workflow_loop = WorkflowLoop(
        rl_agent=rl_agent,
        generative_agent=generative_agent,
        evaluator=evaluator,
        dataset_utils=dataset_utils,
    )
    
    # Initialize Phase 10 orchestrator
    _phase10_orchestrator = Phase10Orchestrator(
        rl_agent=rl_agent,
        generative_agent=generative_agent,
        workflow_loop=workflow_loop,
        orchestrator=_orchestrator,
    )
    
    logger.info("Phase 10 initialized")


# Request/Response models
class GenerateRequest(BaseModel):
    n: int = 10
    method: str = "rl"  # "rl" or "generative"
    seed_smiles: Optional[List[str]] = None


class GenerateResponse(BaseModel):
    molecules: List[str]
    method: str
    count: int


class EvaluateRequest(BaseModel):
    smiles: str
    compute_ml: bool = True
    compute_screening: bool = True
    compute_qm: bool = False
    compute_md: bool = False


class EvaluateResponse(BaseModel):
    smiles: str
    reward: float
    ml_predictions: Dict[str, float]
    screening_results: Dict[str, Any]
    qm_results: Dict[str, float]
    md_results: Dict[str, float]


class RunLoopRequest(BaseModel):
    max_iterations: int = 10
    batch_size: int = 32
    use_generative: bool = False
    seed_smiles: Optional[List[str]] = None


class RunLoopResponse(BaseModel):
    iterations_completed: int
    top_candidates: List[Dict[str, Any]]
    statistics: Dict[str, Any]
    iteration_logs: List[Dict[str, Any]]


class TopCandidatesResponse(BaseModel):
    candidates: List[Dict[str, Any]]
    count: int


@router.on_event("startup")
async def startup_event():
    """Initialize Phase 10 on startup."""
    if ORCHESTRATOR_AVAILABLE:
        try:
            _orchestrator_instance = Orchestrator()
            init_phase10(_orchestrator_instance)
        except Exception as e:
            logger.warning(f"Could not initialize with orchestrator: {e}")
            init_phase10(None)
    else:
        init_phase10(None)


@router.post("/generate", response_model=GenerateResponse)
async def generate_molecules(request: GenerateRequest):
    """
    Generate molecules using RL or generative agent.
    
    Args:
        request: Generation request
    
    Returns:
        Generated molecules
    """
    if _phase10_orchestrator is None:
        raise HTTPException(status_code=503, detail="Phase 10 not initialized")
    
    try:
        molecules = await _phase10_orchestrator.generate_molecules(
            n=request.n,
            method=request.method,
            seed_smiles=request.seed_smiles,
        )
        return GenerateResponse(
            molecules=molecules,
            method=request.method,
            count=len(molecules),
        )
    except Exception as e:
        logger.error(f"Generation error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/evaluate", response_model=EvaluateResponse)
async def evaluate_molecule(request: EvaluateRequest):
    """
    Evaluate a single molecule.
    
    Args:
        request: Evaluation request
    
    Returns:
        Evaluation results with reward
    """
    if _phase10_orchestrator is None:
        raise HTTPException(status_code=503, detail="Phase 10 not initialized")
    
    try:
        evaluator = _phase10_orchestrator.workflow_loop.evaluator
        results = await evaluator.evaluate_molecule(
            smiles=request.smiles,
            compute_ml=request.compute_ml,
            compute_screening=request.compute_screening,
            compute_qm=request.compute_qm,
            compute_md=request.compute_md,
        )
        return EvaluateResponse(
            smiles=results["smiles"],
            reward=results["reward"],
            ml_predictions=results.get("ml_predictions", {}),
            screening_results=results.get("screening_results", {}),
            qm_results=results.get("qm_results", {}),
            md_results=results.get("md_results", {}),
        )
    except Exception as e:
        logger.error(f"Evaluation error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/run_loop", response_model=RunLoopResponse)
async def run_workflow_loop(request: RunLoopRequest):
    """
    Run the RL workflow loop.
    
    Args:
        request: Workflow loop request
    
    Returns:
        Workflow results with top candidates
    """
    if _phase10_orchestrator is None:
        raise HTTPException(status_code=503, detail="Phase 10 not initialized")
    
    try:
        results = await _phase10_orchestrator.run_workflow(
            max_iterations=request.max_iterations,
            batch_size=request.batch_size,
            use_generative=request.use_generative,
            seed_smiles=request.seed_smiles,
        )
        return RunLoopResponse(**results)
    except Exception as e:
        logger.error(f"Workflow loop error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/top_candidates", response_model=TopCandidatesResponse)
async def get_top_candidates(n: int = 10):
    """
    Get top candidates by reward.
    
    Args:
        n: Number of top candidates
    
    Returns:
        Top candidates
    """
    if _phase10_orchestrator is None:
        raise HTTPException(status_code=503, detail="Phase 10 not initialized")
    
    try:
        candidates = _phase10_orchestrator.get_top_candidates(n=n)
        return TopCandidatesResponse(
            candidates=candidates,
            count=len(candidates),
        )
    except Exception as e:
        logger.error(f"Get top candidates error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/statistics")
async def get_statistics():
    """
    Get workflow statistics.
    
    Returns:
        Statistics dictionary
    """
    if _phase10_orchestrator is None:
        raise HTTPException(status_code=503, detail="Phase 10 not initialized")
    
    try:
        stats = _phase10_orchestrator.get_statistics()
        return stats
    except Exception as e:
        logger.error(f"Get statistics error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/iteration_logs")
async def get_iteration_logs():
    """
    Get iteration logs.
    
    Returns:
        List of iteration logs
    """
    if _phase10_orchestrator is None:
        raise HTTPException(status_code=503, detail="Phase 10 not initialized")
    
    try:
        logs = _phase10_orchestrator.get_iteration_logs()
        return {"logs": logs, "count": len(logs)}
    except Exception as e:
        logger.error(f"Get iteration logs error: {e}")
        raise HTTPException(status_code=500, detail=str(e))

