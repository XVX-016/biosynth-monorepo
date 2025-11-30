FROM python:3.11-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    cmake \
    git \
    libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh
ENV PATH="/opt/conda/bin:$PATH"

# Accept Conda Terms of Service for main and r channels
RUN conda config --set always_yes yes --set changeps1 no && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Install RDKit + Python 3.11 from conda-forge
RUN conda install -y -c conda-forge rdkit=2023.09 python=3.11 && \
    conda clean -afy

# Copy requirements first (needed for pip install)
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

# Copy app code into backend/ subdirectory to match import structure
# Since imports use "from backend.config", we need /app/backend/ structure
COPY . /app/backend/

# Set PYTHONPATH to /app so backend.* imports work
ENV PYTHONPATH=/app
ENV PORT=8000

EXPOSE $PORT
# Use PORT env var with default fallback
CMD ["sh", "-c", "uvicorn backend.app:app --host 0.0.0.0 --port ${PORT:-8000}"]
