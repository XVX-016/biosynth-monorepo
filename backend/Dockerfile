FROM python:3.11-slim

WORKDIR /app

# Install system dependencies for RDKit
RUN apt-get update && apt-get install -y \
    build-essential \
    libboost-dev \
    libboost-system-dev \
    libboost-thread-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements first for better caching
COPY backend/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Copy backend code
COPY backend/ ./backend/
# Copy packages if needed (for engine)
COPY packages/ ./packages/

# Set working directory to backend
WORKDIR /app/backend

# Railway provides PORT environment variable
# Default to 8000 for local development
ENV PORT=8000

# Add parent directory to Python path so imports work
ENV PYTHONPATH=/app

CMD uvicorn app:app --host 0.0.0.0 --port ${PORT:-8000}

