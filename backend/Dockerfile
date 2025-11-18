# Optimized Dockerfile using official RDKit base image
# Reduces image size from 6.4GB → ~700MB
# Build context MUST be set to /backend in Railway Settings → Deployment → Directory

# Use official RDKit image with Python 3.11 (latest, includes RDKit + Boost pre-compiled)
FROM rdkit/rdkit:latest-py3.11

WORKDIR /app

# Copy requirements first for better Docker layer caching
COPY requirements.txt .

# Upgrade pip and install Python dependencies
RUN pip install --no-cache-dir --upgrade pip setuptools wheel && \
    pip install --no-cache-dir -r requirements.txt

# Copy backend code (context is already /backend, so no backend/ prefix needed)
COPY . .

# Set Python path for imports
ENV PYTHONPATH=/app
ENV PYTHONUNBUFFERED=1

# Railway provides PORT environment variable
ENV PORT=8000

# Expose port
EXPOSE ${PORT}

# Health check (you have /health endpoint, so this is safe)
HEALTHCHECK --interval=30s --timeout=10s --start-period=40s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:${PORT:-8000}/health')" || exit 1

# Run the application
CMD ["sh", "-c", "uvicorn app:app --host 0.0.0.0 --port ${PORT:-8000}"]
