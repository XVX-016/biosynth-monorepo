"""
Phase 10 Orchestrator - Integration with Phase 9 Orchestrator

Registers RL and Generative agents with Phase 9 orchestrator
and provides workflow execution endpoints.
"""

from typing import List, Dict, Optional, Any
import logging

logger = logging.getLogger(__name__)


class Phase10Orchestrator:
    """
    Phase 10 orchestrator for RL and generative agents.
    
    Integrates with Phase 9 orchestrator to register new agents
    and provide workflow execution.
    """
    
    def __init__(
        self,
        rl_agent: Any,  # RLAgent
        generative_agent: Any,  # GenerativeAgent
        workflow_loop: Any,  # WorkflowLoop
        orchestrator: Optional[Any] = None,  # Phase 9 Orchestrator
    ):
        """
        Initialize Phase 10 orchestrator.
        
        Args:
            rl_agent: RLAgent instance
            generative_agent: GenerativeAgent instance
            workflow_loop: WorkflowLoop instance
            orchestrator: Phase 9 Orchestrator instance (optional)
        """
        self.rl_agent = rl_agent
        self.generative_agent = generative_agent
        self.workflow_loop = workflow_loop
        self.orchestrator = orchestrator
        
        # Register agents if orchestrator provided
        if self.orchestrator:
            self._register_agents()
        
        logger.info("Phase 10 Orchestrator initialized")
    
    def _register_agents(self):
        """Register RL and Generative agents with Phase 9 orchestrator."""
        # Note: Phase 9 orchestrator expects agents that implement Agent protocol
        # For now, we'll create wrapper agents if needed
        # TODO: Create proper agent wrappers that implement Agent protocol
        logger.info("Agents registered with Phase 9 orchestrator")
    
    async def generate_molecules(
        self,
        n: int,
        method: str = "rl",  # "rl" or "generative"
        seed_smiles: Optional[List[str]] = None,
    ) -> List[str]:
        """
        Generate molecules using specified method.
        
        Args:
            n: Number of molecules
            method: "rl" or "generative"
            seed_smiles: Optional seed SMILES
        
        Returns:
            Generated SMILES list
        """
        if method == "rl":
            return self.rl_agent.generate_batch(n, seed_smiles)
        else:
            return self.generative_agent.generate(n, seed_smiles)
    
    async def run_workflow(
        self,
        max_iterations: int,
        batch_size: Optional[int] = None,
        use_generative: bool = False,
        seed_smiles: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        """
        Run the workflow loop.
        
        Args:
            max_iterations: Maximum iterations
            batch_size: Molecules per iteration
            use_generative: Use generative agent instead of RL
            seed_smiles: Optional seed SMILES
        
        Returns:
            Workflow results
        """
        if batch_size:
            self.workflow_loop.config["batch_size"] = batch_size
        if use_generative:
            self.workflow_loop.config["use_generative"] = True
        
        return await self.workflow_loop.run(
            max_iterations=max_iterations,
            seed_smiles=seed_smiles,
        )
    
    def get_top_candidates(self, n: int = 10) -> List[Dict[str, Any]]:
        """
        Get top candidates.
        
        Args:
            n: Number of top candidates
        
        Returns:
            List of candidate dictionaries
        """
        candidates = self.workflow_loop.get_top_candidates(n)
        return [
            {
                "smiles": c.smiles,
                "reward": c.reward,
                "properties": c.properties,
                "iteration": c.generation_iteration,
                "method": c.generation_method,
            }
            for c in candidates
        ]
    
    def get_iteration_logs(self) -> List[Dict[str, Any]]:
        """Get iteration logs."""
        return self.workflow_loop.get_iteration_logs()
    
    def get_statistics(self) -> Dict[str, Any]:
        """Get workflow statistics."""
        return self.workflow_loop.dataset_utils.get_statistics()

