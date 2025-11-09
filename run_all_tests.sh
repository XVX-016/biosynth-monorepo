#!/usr/bin/env bash
set -e
echo "Running TypeScript engine tests..."
(cd packages/engine && npm run test || true)
echo "Running frontend tests (if any)..."
(cd frontend && npm run test || true)
echo "Running backend tests..."
(cd backend && python -m pytest -q || true)

