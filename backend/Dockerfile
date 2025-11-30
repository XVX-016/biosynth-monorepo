# === Base Image ===
FROM python:3.11-slim

# === Set Working Directory ===
WORKDIR /app

# === Install System Dependencies ===
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    cmake \
    git \
    libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*

# === Install Miniconda ===
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh
ENV PATH="/opt/conda/bin:$PATH"

# === Conda Config & Accept Terms ===
RUN conda config --set always_yes yes --set changeps1 no && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# === Install RDKit + Python 3.11 ===
RUN conda install -y -c conda-forge rdkit=2023.09 python=3.11 && \
    conda clean -afy

# === Copy & Install Python Dependencies ===
# Remove rdkit-pypi from requirements.txt before this step
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

# === Copy App Code ===
COPY . /app/backend/

# === Environment Variables ===
ENV PYTHONPATH=/app
ENV PORT=8000

# === Expose Port ===
EXPOSE $PORT

# === Entrypoint ===
CMD ["sh", "-c", "uvicorn backend.app:app --host 0.0.0.0 --port ${PORT:-8000}"]
