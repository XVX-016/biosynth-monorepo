FROM rdkit/rdkit:2023.09.1-py3.11-slim

WORKDIR /app

# Upgrade pip, install CPU PyTorch
COPY requirements-model.txt .
RUN pip install --no-cache-dir --upgrade pip setuptools wheel && \
    pip install --no-cache-dir torch==2.4.1+cpu \
        torchvision==0.19.1+cpu \
        torchaudio==2.4.1+cpu \
        --index-url https://download.pytorch.org/whl/cpu && \
    pip install --no-cache-dir -r requirements-model.txt

# Install lightweight dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy backend code (context is already /backend, so no backend/ prefix needed)
COPY . .

ENV PYTHONPATH=/app
ENV PYTHONUNBUFFERED=1
# PORT is provided by platform (Railway or Cloud Run)
# Railway uses $PORT, Cloud Run uses 8080
ENV PORT=8080

EXPOSE ${PORT}

HEALTHCHECK --interval=30s --timeout=10s --start-period=40s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:${PORT:-8080}/health')" || exit 1

CMD ["sh", "-c", "uvicorn app:app --host 0.0.0.0 --port ${PORT:-8080}"]
