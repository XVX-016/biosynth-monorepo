# Phase 5 ML Model Training Pipeline

## Overview

Complete training pipeline for Phase 5 Attention-GNN models using prepared datasets.

## Pipeline Steps

### Step 1: Load Datasets (`load_datasets.py`)
- Loads training, validation, and test datasets from `/data/datasets/cleaned/`
- Loads associated fingerprints from `/data/datasets/fingerprints/`
- Ensures SMILES and property values align correctly

### Step 2: Configure Models (`configure_model.py`)
- Initializes Attention-GNN models for each property (logP, solubility, toxicity)
- Sets input dimensions from sample molecule featurization
- Configures hyperparameters (hidden_dim, num_layers, heads, etc.)

### Step 3: Train Models (`train_models.py`)
- Trains models using training set (70%)
- Validates on validation set (15%)
- Tracks loss curves and saves best checkpoints
- Saves models to `/data/models/{property}_best.pt`

### Step 4: Evaluate Models (`evaluate_models.py`)
- Runs predictions on test set (15%)
- Computes metrics: MAE, RMSE, R²
- Saves evaluation reports to `/data/models/evaluation/`

### Step 5: Register Models (`register_models.py`)
- Registers trained models in Phase 5 model registry
- Includes metadata: training date, dataset, properties, hyperparameters, metrics
- Ensures API endpoints can load and use these models

### Step 6: Test API (`test_predictions.py`)
- Tests `POST /api/predict/property` endpoint
- Confirms output values are realistic
- Confirms attention maps are generated (if applicable)

### Step 7: Precompute Predictions (`precompute_predictions.py`)
- Precomputes property predictions for library molecules
- Reads from `/data/libraries/*.smi`
- Stores results in `/data/libraries/predictions.csv`
- Used by Phase 10 RL loop and screening agents

## Usage

### Run Full Pipeline

```bash
cd backend
python scripts/train_ml_models.py
```

### Run Individual Steps

```bash
# Step 1: Load datasets
python scripts/load_datasets.py

# Step 2: Configure models
python scripts/configure_model.py

# Step 3: Train models
python scripts/train_models.py

# Step 4: Evaluate models
python scripts/evaluate_models.py

# Step 5: Register models
python scripts/register_models.py

# Step 6: Test API
python scripts/test_predictions.py

# Step 7: Precompute predictions
python scripts/precompute_predictions.py
```

## Requirements

- PyTorch
- PyTorch Geometric
- RDKit
- NumPy
- FastAPI (for API testing)

## Model Architecture

- **Type**: Attention-GNN (Graph Attention Network)
- **Input**: Node features (64D), Edge features (8D)
- **Architecture**: 3 GAT layers, 4 attention heads, 128 hidden dim
- **Output**: Single property value per model

## Output Files

- **Models**: `/data/models/{property}_best.pt`
- **Registry**: `/data/models/registry.json`
- **Evaluations**: `/data/models/evaluation/{property}_evaluation.json`
- **Precomputed**: `/data/libraries/predictions.csv`

## Notes

- Models are trained separately for each property
- Training uses MSE loss and Adam optimizer
- Best model is saved based on validation loss
- Models are registered with metadata for API access
- Precomputed predictions speed up Phase 10 RL loop

## Integration

After training:
1. Models are automatically registered in ModelRegistry
2. API endpoints (`/api/predict/property`) can use trained models
3. Phase 10 RL loop can use real predictions instead of mocks
4. Screening agents can use precomputed predictions

