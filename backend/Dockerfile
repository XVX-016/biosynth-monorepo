FROM python:3.10-slim

WORKDIR /app

# Install system dependencies for RDKit
RUN apt-get update && apt-get install -y \
    build-essential \
    libboost-dev \
    libboost-system-dev \
    libboost-thread-dev \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "8000"]

