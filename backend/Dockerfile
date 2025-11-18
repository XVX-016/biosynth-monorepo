FROM rdkit/rdkit:2023.09.1-py3.11-slim

WORKDIR /app

# Upgrade pip first
RUN pip install --no-cache-dir --upgrade pip setuptools wheel

# Install CPU-only PyTorch FIRST (before other deps to avoid conflicts)
# Using CPU-only wheels reduces size from ~800MB to ~120MB
# This MUST be installed separately to ensure CPU-only version
RUN pip install --no-cache-dir \
    torch==2.4.1+cpu \
    torchvision==0.19.1+cpu \
    torchaudio==2.4.1+cpu \
    --index-url https://download.pytorch.org/whl/cpu

# Install ML dependencies (transformers, onnxruntime, etc.)
# Build context is repo root, so use backend/ prefix
COPY backend/requirements-model.txt .
RUN pip install --no-cache-dir -r requirements-model.txt

# Install lightweight dependencies
COPY backend/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy backend code (build context is repo root, so use backend/ prefix)
COPY backend/ .

ENV PYTHONPATH=/app
ENV PYTHONUNBUFFERED=1
# PORT is provided by platform (Railway or Cloud Run)
# Railway uses $PORT, Cloud Run uses 8080
ENV PORT=8080

EXPOSE ${PORT}

HEALTHCHECK --interval=30s --timeout=10s --start-period=40s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:${PORT:-8080}/health')" || exit 1

CMD ["sh", "-c", "uvicorn app:app --host 0.0.0.0 --port ${PORT:-8080}"]
