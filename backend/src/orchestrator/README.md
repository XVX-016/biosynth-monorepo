# Phase 9: Multi-Agent Drug Discovery Orchestrator

## Overview

Phase 9 implements a multi-agent orchestration system for coordinating drug discovery workflows. Agents can execute tasks independently and workflows can chain multiple tasks together.

## Structure

```
backend/src/orchestrator/
├── __init__.py              # Module exports
├── agent_protocols.py       # Base Agent class, Task, TaskResult, Message
├── agent_registry.py         # AgentRegistry for managing agents
├── task_router.py           # TaskRouter for routing tasks to agents
├── orchestrator.py           # Main Orchestrator class
├── workflow_specs.json      # Predefined workflow templates
└── agents/
    ├── __init__.py
    ├── predictor_agent.py   # ML prediction agent
    ├── screening_agent.py   # Screening agent
    ├── qm_agent.py          # Quantum chemistry agent
    └── md_agent.py          # Molecular dynamics agent
```

## Components

### 1. Agent Protocols (`agent_protocols.py`)

**Base Classes:**
- `Agent` - Abstract base class for all agents
- `AgentProtocol` - Protocol interface
- `Task` - Task definition
- `TaskResult` - Task execution result
- `Message` - JSON message protocol

**Key Methods:**
- `Agent.get_name()` - Get agent name
- `Agent.get_config()` - Get agent configuration
- `Agent.run(task)` - Execute task (async)

### 2. Agent Registry (`agent_registry.py`)

**Features:**
- Register/unregister agents
- Track agent capabilities
- Route tasks to capable agents
- Get statistics

**Usage:**
```python
registry = AgentRegistry()
registry.register(PredictorAgent())
agents = registry.get_agents_for_task('predict')
```

### 3. Task Router (`task_router.py`)

**Routing Strategies:**
- `round_robin` - Distribute tasks evenly
- `rule_based` - Route based on task properties
- `random` - Random selection

**Usage:**
```python
router = TaskRouter(registry)
agent = router.route(task, strategy='round_robin')
```

### 4. Orchestrator (`orchestrator.py`)

**Features:**
- Submit tasks
- Execute tasks
- Execute workflows (sequence of tasks)
- Track task status and history
- Get statistics

**Usage:**
```python
orchestrator = Orchestrator()
orchestrator.register_agent(PredictorAgent())

task_id = orchestrator.submit_task(
    task_type='predict',
    input_data={'smiles': 'CCO'}
)
result = await orchestrator.execute_task(task_id)
```

### 5. Mock Agents

#### Predictor Agent (`agents/predictor_agent.py`)
- Integrates with Phase 5 ML engine
- Handles `predict` tasks
- Returns property predictions

#### Screening Agent (`agents/screening_agent.py`)
- Integrates with Phase 7 search engine
- Handles `screen` tasks
- Supports similarity and substructure search

#### QM Agent (`agents/qm_agent.py`)
- Integrates with Phase 8 QM engine
- Handles `qm_energy`, `qm_optimize`, `qm_properties` tasks
- Uses mock QM calculations

#### MD Agent (`agents/md_agent.py`)
- Integrates with Phase 8 MD engine
- Handles `md_simulate` tasks
- Uses mock MD simulations

## Message Protocol

All communication uses JSON messages:

```json
{
  "message_id": "msg_123",
  "message_type": "task",
  "sender": "orchestrator",
  "receiver": "predictor",
  "payload": {
    "task_id": "task_456",
    "input_data": {"smiles": "CCO"}
  },
  "timestamp": 1234567890.0,
  "metadata": {}
}
```

## API Endpoints

### Task Management
- `POST /api/orchestrator/task/submit` - Submit a task
- `POST /api/orchestrator/task/execute` - Execute a task
- `GET /api/orchestrator/task/status/{task_id}` - Get task status

### Workflow Execution
- `POST /api/orchestrator/workflow/execute` - Execute a workflow

### Agent Management
- `GET /api/orchestrator/agents` - List all agents
- `GET /api/orchestrator/stats` - Get orchestrator statistics

## Usage Examples

### Submit and Execute Task

```python
POST /api/orchestrator/task/submit
{
  "task_type": "predict",
  "input_data": {
    "smiles": "CCO",
    "properties": ["logP", "toxicity"]
  },
  "priority": 1
}

# Returns: {"task_id": "..."}

POST /api/orchestrator/task/execute?task_id=...
# Returns: TaskResult with predictions
```

### Execute Workflow

```python
POST /api/orchestrator/workflow/execute
{
  "workflow": [
    {
      "task_type": "screen",
      "input_data": {
        "task_subtype": "similarity",
        "query_smiles": "CCO",
        "k": 10
      },
      "priority": 2
    },
    {
      "task_type": "predict",
      "input_data": {
        "smiles": "CCO",
        "properties": ["logP"]
      },
      "priority": 1
    }
  ],
  "routing_strategy": "round_robin"
}
```

## Workflow Specs

Predefined workflows are in `workflow_specs.json`:
- `predict_properties` - Simple property prediction
- `screen_and_predict` - Screen library then predict
- `optimize_and_simulate` - QM optimize then MD simulate
- `full_pipeline` - Complete discovery pipeline

## Future Enhancements

TODO comments left for:
- Real docking agent integration
- Real MD engine integration (not mock)
- Generative model agents
- Retrosynthesis agents
- Load-based routing
- Task prioritization queue
- Parallel task execution
- Workflow template engine

## Integration

- **Phase 5**: PredictorAgent uses ML engine
- **Phase 7**: ScreeningAgent uses search engine
- **Phase 8**: QMAgent and MDAgent use QM/MD engines

All agents integrate via Python imports (no external calls).

