"""
Integration tests for Phase 10: RL + Generative Molecule Design

Tests:
- RL agent generation
- Generative agent generation
- Reward function computation
- Evaluator integration
- Workflow loop execution
- Top candidates endpoint
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

from src.phase10 import (
    RLAgent,
    GenerativeAgent,
    RewardFunction,
    Evaluator,
    WorkflowLoop,
    DatasetUtils,
    Phase10Orchestrator,
)


class TestPhase10RLAgent:
    """Test RL Agent."""
    
    def test_rl_agent_init(self):
        """Test RL agent initialization."""
        agent = RLAgent()
        assert agent is not None
        assert agent.policy_state.generation_count == 0
    
    def test_rl_agent_generate_batch(self):
        """Test batch generation."""
        agent = RLAgent()
        molecules = agent.generate_batch(n=5)
        assert len(molecules) == 5
        assert all(isinstance(m, str) for m in molecules)
        assert agent.policy_state.generation_count == 5
    
    def test_rl_agent_update_policy(self):
        """Test policy update."""
        agent = RLAgent()
        molecules = ["CCO", "CCCO", "CC(C)O"]
        rewards = [0.5, 0.7, 0.3]
        
        update = agent.update_policy(rewards, molecules)
        assert update["updated"] is True
        assert update["max_reward"] == 0.7
        assert agent.policy_state.best_reward == 0.7


class TestPhase10GenerativeAgent:
    """Test Generative Agent."""
    
    def test_generative_agent_init(self):
        """Test generative agent initialization."""
        agent = GenerativeAgent()
        assert agent is not None
        assert agent.config["model_type"] == "diffusion"
    
    def test_generative_agent_generate(self):
        """Test molecule generation."""
        agent = GenerativeAgent()
        molecules = agent.generate(n=5)
        assert len(molecules) == 5
        assert all(isinstance(m, str) for m in molecules)
    
    def test_generative_agent_with_seeds(self):
        """Test generation with seed SMILES."""
        agent = GenerativeAgent()
        seeds = ["CCO", "CCCO"]
        molecules = agent.generate(n=4, seed_smiles=seeds)
        assert len(molecules) == 4


class TestPhase10RewardFunction:
    """Test Reward Function."""
    
    def test_reward_function_init(self):
        """Test reward function initialization."""
        rf = RewardFunction()
        assert rf is not None
        assert "logP" in rf.config["weights"]
    
    def test_reward_function_compute(self):
        """Test reward computation."""
        rf = RewardFunction()
        ml_predictions = {"logP": 2.5, "solubility": 0.8, "toxicity": 0.1}
        reward = rf.compute("CCO", ml_predictions=ml_predictions)
        assert isinstance(reward, float)
    
    def test_reward_function_batch_compute(self):
        """Test batch reward computation."""
        rf = RewardFunction()
        smiles_list = ["CCO", "CCCO", "CC(C)O"]
        predictions_list = [
            {"logP": 2.5},
            {"logP": 3.0},
            {"logP": 2.0},
        ]
        rewards = rf.batch_compute(smiles_list, predictions_list=predictions_list)
        assert len(rewards) == 3
        assert all(isinstance(r, float) for r in rewards)


class TestPhase10DatasetUtils:
    """Test Dataset Utils."""
    
    def test_dataset_utils_init(self):
        """Test dataset utils initialization."""
        utils = DatasetUtils()
        assert utils is not None
    
    def test_add_record(self):
        """Test adding records."""
        utils = DatasetUtils()
        record = utils.add_record(
            smiles="CCO",
            reward=0.5,
            properties={"logP": 2.5},
            generation_iteration=1,
            generation_method="rl",
        )
        assert record.smiles == "CCO"
        assert len(utils.records) == 1
    
    def test_get_top_candidates(self):
        """Test getting top candidates."""
        utils = DatasetUtils()
        utils.add_record("CCO", 0.3, {}, 1, "rl")
        utils.add_record("CCCO", 0.7, {}, 1, "rl")
        utils.add_record("CC(C)O", 0.5, {}, 1, "rl")
        
        top = utils.get_top_candidates(n=2)
        assert len(top) == 2
        assert top[0].reward == 0.7
        assert top[1].reward == 0.5


class TestPhase10Evaluator:
    """Test Evaluator."""
    
    @pytest.mark.asyncio
    async def test_evaluator_init(self):
        """Test evaluator initialization."""
        rf = RewardFunction()
        evaluator = Evaluator(reward_function=rf)
        assert evaluator is not None
    
    @pytest.mark.asyncio
    async def test_evaluate_molecule(self):
        """Test molecule evaluation."""
        rf = RewardFunction()
        evaluator = Evaluator(reward_function=rf)
        
        # Mock evaluation (without orchestrator)
        results = await evaluator.evaluate_molecule(
            "CCO",
            compute_ml=False,  # Skip ML to avoid orchestrator dependency
            compute_screening=False,
        )
        assert results["smiles"] == "CCO"
        assert "reward" in results


class TestPhase10WorkflowLoop:
    """Test Workflow Loop."""
    
    @pytest.mark.asyncio
    async def test_workflow_loop_init(self):
        """Test workflow loop initialization."""
        rl_agent = RLAgent()
        gen_agent = GenerativeAgent()
        rf = RewardFunction()
        evaluator = Evaluator(reward_function=rf)
        dataset_utils = DatasetUtils()
        
        loop = WorkflowLoop(
            rl_agent=rl_agent,
            generative_agent=gen_agent,
            evaluator=evaluator,
            dataset_utils=dataset_utils,
        )
        assert loop is not None
    
    @pytest.mark.asyncio
    async def test_run_iteration(self):
        """Test running a single iteration."""
        rl_agent = RLAgent()
        gen_agent = GenerativeAgent()
        rf = RewardFunction()
        evaluator = Evaluator(reward_function=rf)
        dataset_utils = DatasetUtils()
        
        loop = WorkflowLoop(
            rl_agent=rl_agent,
            generative_agent=gen_agent,
            evaluator=evaluator,
            dataset_utils=dataset_utils,
            config={"batch_size": 5, "max_iterations": 10},
        )
        
        result = await loop.run_iteration(1)
        assert result["iteration"] == 1
        assert "avg_reward" in result


class TestPhase10Orchestrator:
    """Test Phase 10 Orchestrator."""
    
    @pytest.mark.asyncio
    async def test_phase10_orchestrator_init(self):
        """Test Phase 10 orchestrator initialization."""
        rl_agent = RLAgent()
        gen_agent = GenerativeAgent()
        rf = RewardFunction()
        evaluator = Evaluator(reward_function=rf)
        dataset_utils = DatasetUtils()
        loop = WorkflowLoop(
            rl_agent=rl_agent,
            generative_agent=gen_agent,
            evaluator=evaluator,
            dataset_utils=dataset_utils,
        )
        
        orchestrator = Phase10Orchestrator(
            rl_agent=rl_agent,
            generative_agent=gen_agent,
            workflow_loop=loop,
        )
        assert orchestrator is not None
    
    @pytest.mark.asyncio
    async def test_generate_molecules(self):
        """Test molecule generation via orchestrator."""
        rl_agent = RLAgent()
        gen_agent = GenerativeAgent()
        rf = RewardFunction()
        evaluator = Evaluator(reward_function=rf)
        dataset_utils = DatasetUtils()
        loop = WorkflowLoop(
            rl_agent=rl_agent,
            generative_agent=gen_agent,
            evaluator=evaluator,
            dataset_utils=dataset_utils,
        )
        
        orchestrator = Phase10Orchestrator(
            rl_agent=rl_agent,
            generative_agent=gen_agent,
            workflow_loop=loop,
        )
        
        molecules = await orchestrator.generate_molecules(n=5, method="rl")
        assert len(molecules) == 5


class TestPhase10API:
    """Test API imports."""
    
    def test_api_imports(self):
        """Test that API module can be imported."""
        try:
            from backend.api import phase10
            assert phase10.router is not None
        except ImportError as e:
            pytest.skip(f"API module not available: {e}")

