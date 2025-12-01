# Phase 10: RL + Generative Molecule Design

## Overview

Phase 10 implements reinforcement learning (RL) and generative model-based molecule design. The system generates molecules, evaluates them using ML predictions, screening, and QM/MD calculations, computes rewards, and updates generation policies iteratively.

## Architecture

### Components

1. **RLAgent** (`rl_agent.py`)
   - Policy-based molecule generation
   - Batch generation with seed molecules
   - Policy updates based on reward feedback
   - Tracks generation history and best rewards

2. **GenerativeAgent** (`generative_agent.py`)
   - Diffusion/VAE-based molecule generation
   - Supports seed-based and unconditional generation
   - Configurable model types and temperature

3. **RewardFunction** (`reward_function.py`)
   - Composite reward scoring
   - Combines ML predictions, screening, QM/MD results
   - Configurable property weights
   - Constraint application

4. **Evaluator** (`evaluator.py`)
   - Orchestrates molecule evaluation
   - Integrates with Phase 5 (ML), Phase 7 (screening), Phase 8 (QM/MD)
   - Batch evaluation support
   - Returns comprehensive evaluation results

5. **WorkflowLoop** (`workflow_loop.py`)
   - Main optimization loop
   - Batch generation → evaluation → policy update
   - Iteration logging and top candidate tracking
   - Early stopping support

6. **DatasetUtils** (`dataset_utils.py`)
   - Manages generated molecule datasets
   - Stores records with rewards and properties
   - Top candidate retrieval
   - Dataset persistence (JSON)

7. **Phase10Orchestrator** (`phase10_orchestrator.py`)
   - High-level interface for Phase 10
   - Integrates with Phase 9 orchestrator
   - Workflow execution management

## API Endpoints

### POST `/api/phase10/generate`
Generate molecules using RL or generative agent.

**Request:**
```json
{
  "n": 10,
  "method": "rl",  // or "generative"
  "seed_smiles": ["CCO", "CCCO"]  // optional
}
```

**Response:**
```json
{
  "molecules": ["CCO", "CCCO", ...],
  "method": "rl",
  "count": 10
}
```

### POST `/api/phase10/evaluate`
Evaluate a single molecule.

**Request:**
```json
{
  "smiles": "CCO",
  "compute_ml": true,
  "compute_screening": true,
  "compute_qm": false,
  "compute_md": false
}
```

**Response:**
```json
{
  "smiles": "CCO",
  "reward": 0.75,
  "ml_predictions": {"logP": 2.5, "solubility": 0.8},
  "screening_results": {},
  "qm_results": {},
  "md_results": {}
}
```

### POST `/api/phase10/run_loop`
Run the RL workflow loop.

**Request:**
```json
{
  "max_iterations": 10,
  "batch_size": 32,
  "use_generative": false,
  "seed_smiles": ["CCO"]  // optional
}
```

**Response:**
```json
{
  "iterations_completed": 10,
  "top_candidates": [
    {
      "smiles": "CCO",
      "reward": 0.85,
      "properties": {"logP": 2.5},
      "iteration": 5
    }
  ],
  "statistics": {...},
  "iteration_logs": [...]
}
```

### GET `/api/phase10/top_candidates?n=10`
Get top candidates by reward.

### GET `/api/phase10/statistics`
Get workflow statistics.

### GET `/api/phase10/iteration_logs`
Get iteration logs.

## Integration

### Phase 5 (ML Engine)
- Uses ML predictions for reward computation
- Properties: logP, solubility, toxicity, etc.

### Phase 7 (Search & Screening)
- Uses similarity search for reward computation
- Screening results contribute to reward

### Phase 8 (QM/MD)
- Optional QM energy calculations
- Optional MD simulations
- Results contribute to reward

### Phase 9 (Orchestrator)
- Evaluator uses orchestrator to execute tasks
- Agents can be registered with orchestrator

## Usage Example

```python
from src.phase10 import (
    RLAgent,
    GenerativeAgent,
    RewardFunction,
    Evaluator,
    WorkflowLoop,
    DatasetUtils,
)

# Initialize components
rl_agent = RLAgent()
gen_agent = GenerativeAgent()
rf = RewardFunction()
evaluator = Evaluator(reward_function=rf)
dataset_utils = DatasetUtils()

# Create workflow loop
loop = WorkflowLoop(
    rl_agent=rl_agent,
    generative_agent=gen_agent,
    evaluator=evaluator,
    dataset_utils=dataset_utils,
)

# Run optimization
results = await loop.run(max_iterations=10)

# Get top candidates
top = loop.get_top_candidates(n=5)
```

## Configuration

### RL Agent
- `learning_rate`: Policy update learning rate (default: 0.001)
- `batch_size`: Molecules per batch (default: 32)
- `exploration_rate`: Exploration vs exploitation (default: 0.1)

### Generative Agent
- `model_type`: "diffusion" or "vae" (default: "diffusion")
- `seed_library`: Path to seed library (optional)
- `temperature`: Sampling temperature (default: 1.0)

### Reward Function
- `weights`: Property weights (logP, solubility, toxicity, etc.)
- `constraints`: Constraint functions
- `normalization`: Whether to normalize scores (default: True)

### Workflow Loop
- `batch_size`: Molecules per iteration (default: 32)
- `max_iterations`: Maximum iterations (default: 100)
- `use_generative`: Use generative agent instead of RL (default: False)
- `top_k`: Number of top candidates to track (default: 10)

## Future Enhancements

1. **Real RL Implementation**
   - Policy network (GNN-based)
   - Value function
   - Actor-critic architecture

2. **Real Generative Models**
   - Diffusion model training
   - VAE training
   - Pre-trained model loading

3. **Advanced Reward Functions**
   - Multi-objective optimization
   - Pareto frontier
   - Constraint satisfaction

4. **Active Learning**
   - Feed top molecules back to ML training
   - Iterative model improvement

5. **GPU Acceleration**
   - Batch generation on GPU
   - Parallel evaluation

## Testing

Run integration tests:
```bash
pytest backend/tests/test_phase10_integration.py -v
```

## Notes

- Current implementation uses mock generation for testing
- Real RL and generative models require training data and model files
- Integration with Phase 9 orchestrator is optional but recommended
- QM/MD calculations are expensive and optional

