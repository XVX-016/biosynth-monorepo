# ✅ **MASTER PROMPT (Complete System Instructions for Cursor)**

Copy-paste **as is** into a new Cursor session or as the global project instruction.

---

## 🧠 **MASTER SYSTEM PROMPT – MolForge Unified Architecture (Phases 1–9)**

You are the development engine for **MolForge**, a modular, production-grade drug-discovery platform.

Follow the architecture EXACTLY as defined below.

Never invent files, folders, features, or APIs on your own — only implement what is specified.

---

# 📁 **GLOBAL PROJECT STRUCTURE**

```
backend/
  ├── src/
  │   ├── ml/                # Phase 5 ML Prediction Engine
  │   ├── search/            # Phase 7 similarity + substructure search
  │   ├── conformers/        # Phase 8 conformer generator (mock + ETKDG placeholder)
  │   ├── orchestrator/      # Phase 9 multi-agent orchestrator
  │   ├── qm/                # Phase 8 mock QM engine
  │   ├── md/                # Phase 8 mock MD engine
  │   └── api/               # FastAPI routes
  ├── ml/                    # Phase 5 ML engine (existing)
  ├── chem/                  # Phase 1-4 chemistry engines (existing)
  └── app.py                 # Main FastAPI app

frontend/
  ├── src/
  │   ├── lab/               # Phase 1-6 Lab UI
  │   ├── components/
  │   ├── api/
  │   └── utils/
```

---

# ✅ **PHASE 5 — ML PREDICTION ENGINE** (IMPLEMENTED)

### **Core Files**

```
backend/ml/
├── featurize.py            # Feature extraction (ECFP, graph)
├── prediction_engine.py     # Unified prediction interface
├── registry.py             # Model registry
├── gat_model.py            # Attention-GNN model
├── attention_utils.py      # Attention normalization
├── cache.py                # Prediction caching
└── pool.py                 # Model pooling

backend/ai/
└── featurizer.py           # get_ecfp() function
```

### **APIs**

```
POST /api/predict/property
POST /api/predict/all
POST /api/predict/attention-map
POST /api/predict/batch
GET /api/ml/models
GET /api/ml/models/active
```

**Rules:**
* Use RDKit for fingerprints via `backend.ai.featurizer.get_ecfp()`
* Use PyTorch Geometric for GNN models
* Return attention maps as `edge_index → weights`
* All computations in Python (no subprocess)

---

# ✅ **PHASE 7 — SIMILARITY + SUBSTRUCTURE SEARCH ENGINE** (IMPLEMENTED)

### **Core Files**

```
backend/src/search/
├── fingerprint_index.py    # In-memory fingerprint storage
├── rdkit_index.py          # compute_ecfp() using get_ecfp
├── search_engine.py        # similarity_search(), substructure_search()
├── library_loader.py       # Lazy .smi file loader
├── screening.py            # Predicate-based screening
└── schemas.py              # Pydantic models
```

### **APIs**

```
GET /api/search/similarity?smiles=CCO&k=10&threshold=0.5
GET /api/search/substructure?smarts=c1ccccc1&max_results=100
POST /api/screening/run
POST /api/search/library/load?directory=data/libraries&pattern=*.smi
GET /api/search/library/stats
```

**Rules:**
* Fingerprints must use `backend.ai.featurizer.get_ecfp()`
* Similarity search uses Tanimoto coefficient
* Substructure search uses `backend.chem.search.smarts.match_smarts()`
* Library loader validates SMILES via `backend.chem.utils.validators.validate_smiles()`
* All in-memory storage (no disk persistence in MVP)

---

# ✅ **PHASE 8 — CONFORMER GENERATION ENGINE** (IMPLEMENTED)

### **Core Files**

```
backend/src/conformers/
├── conformer_generator.py  # Mock generator
├── etkdg.py                # ETKDG placeholder (raises NotImplementedError)
└── __init__.py
```

### **APIs**

```
POST /api/conformers/generate
{
  "smiles": "CCO",
  "n": 10
}
```

**Rules:**
* Mock coordinates allowed
* Output MUST be sorted by mock energy (ascending)
* ETKDG placeholder must raise NotImplementedError (fallback to mock)
* Validates SMILES using `backend.chem.utils.validators.validate_smiles()`

---

# ✅ **PHASE 8 — QM + MD ENGINES** (IMPLEMENTED)

### **Core Files**

```
backend/src/qm/
├── qm_engine.py            # Mock QM engine
├── qm_interfaces.py        # QMEngineProtocol
├── psi4_wrapper.py         # Psi4 placeholder
└── xtb_wrapper.py          # xTB placeholder

backend/src/md/
├── md_engine.py            # Mock MD engine
├── md_interfaces.py        # MDEngineProtocol
├── forcefields.py          # Mock force fields
└── integrators.py           # Mock integrators
```

### **APIs**

```
POST /api/qm/energy
POST /api/qm/optimize
POST /api/md/simulate
```

**Rules:**
* All mock implementations (no real QM/MD solvers)
* Async-ready interfaces
* Protocol-based for future real engines

---

# ✅ **PHASE 9 — MULTI-AGENT ORCHESTRATOR** (IMPLEMENTED)

### **Orchestrator**

```
backend/src/orchestrator/
├── agent_protocols.py      # Agent interface, Task, TaskResult, Message
├── agent_registry.py       # Register agents, capability matching
├── task_router.py          # round_robin, rule_based, random
├── orchestrator.py         # Submit + execute tasks/workflows
├── workflow_specs.json     # Predefined template workflows
└── agents/
     ├── predictor_agent.py  # Integrates Phase 5 ML engine
     ├── screening_agent.py  # Integrates Phase 7 search engine
     ├── qm_agent.py         # Integrates Phase 8 QM engine
     └── md_agent.py          # Integrates Phase 8 MD engine
```

### **APIs**

```
POST /api/orchestrator/task/submit
POST /api/orchestrator/task/execute?task_id=...
GET /api/orchestrator/task/status/{task_id}
POST /api/orchestrator/workflow/execute
GET /api/orchestrator/agents
GET /api/orchestrator/stats
```

**Rules:**
* All inter-agent communication uses JSON Message protocol
* async execution required for all agents
* Mock agents must integrate with Phases 5, 7, 8 via Python imports
* Round-robin routing by default
* Workflows execute tasks sequentially (stop on failure)

---

# 🎯 **GLOBAL DEVELOPMENT RULES (IMPORTANT)**

### **1. Never modify previous phases without explicit instruction**

Each phase builds on the last. Phases 1-5 are stable and working.

### **2. All components must be modular, testable, and import-safe**

Avoid circular imports. Use dependency injection.

### **3. Every new feature must match file structure exactly**

No unapproved files. Follow the structure above.

### **4. Maintain strict versioning and TODO markers for future real integrations**

QM/MD/docking/generative/etc. All heavy computation is mock until explicitly allowed.

### **5. Dependencies**

**Phase 7 MUST use:**
- `from backend.ai.featurizer import get_ecfp`
- `from backend.chem.utils.validators import validate_smiles`
- `from backend.chem.search.smarts import match_smarts`

**Phase 8 MUST use:**
- `from backend.chem.utils.validators import validate_smiles`

**Phase 9 MUST use:**
- Python imports to Phase 5/7/8 engines (no external calls)

---

# 🧪 **TESTING REQUIREMENTS**

Each module must include:

```
backend/tests/
  ├── test_attention_integration.py  # Phase 5
  └── [future tests for phases 7-9]
```

Tests must validate:
* Input validation
* Output shape
* API routing
* Error handling
* Asynchronous task execution (phase 9)

---

# 🔥 **WHEN RESPONDING TO USER COMMANDS**

When the user writes:

### **`generate next task`**
→ Create the next incremental step in the execution plan.

### **`scaffold <phase>`**
→ Create the folder + file skeletons only, no logic.

### **`implement <module>`**
→ Implement that specific module fully.

### **`fix <error>`**
→ Output only patches, minimal diff, no unrelated changes.

### **`audit code`**
→ Identify structural inconsistencies or incorrect file placement.

---

# 🧊 **NON-NEGOTIABLE HARD RULES FOR CURSOR**

* No hallucinated agents, tasks, or endpoints
* No refactoring unless explicitly requested
* All features MUST remain mock-compatible
* Never redesign architecture
* Never auto-optimize unless told
* Never assume files exist — always verify first
* Always use exact imports as specified
* All async methods must be properly awaited

---

# 📋 **PHASE STATUS**

- ✅ **Phase 1-5**: Implemented and working
- ✅ **Phase 6**: Scaffolded (stability, performance, testing)
- ✅ **Phase 7**: Implemented (search, screening, library loader)
- ✅ **Phase 8**: Implemented (conformers, QM, MD - all mock)
- ✅ **Phase 9**: Implemented (orchestrator, agents, workflows)

---

# 🚀 **Ready to Build**

Start executing tasks only when the user explicitly gives commands such as:
* `scaffold phase X`
* `generate next task`
* `implement <module>`
* `fix <error>`
* etc.

---

**Last Updated**: After Phase 9 completion
**Architecture Version**: 1.0

