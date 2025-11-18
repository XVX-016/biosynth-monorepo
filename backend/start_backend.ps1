# PowerShell script to start the backend server
# Usage: .\start_backend.ps1

Write-Host "Starting MolForge Backend Server..." -ForegroundColor Green

# Activate virtual environment
if (Test-Path ".venv\Scripts\Activate.ps1") {
    .\.venv\Scripts\Activate.ps1
    Write-Host "✅ Virtual environment activated" -ForegroundColor Green
} else {
    Write-Host "❌ Virtual environment not found at .venv\Scripts\Activate.ps1" -ForegroundColor Red
    exit 1
}

# Check if uvicorn is installed
python -c "import uvicorn" 2>$null
if ($LASTEXITCODE -ne 0) {
    Write-Host "❌ uvicorn not installed. Installing..." -ForegroundColor Yellow
    pip install uvicorn
}

# Start the server
Write-Host "🚀 Starting server on http://localhost:8000" -ForegroundColor Cyan
Write-Host "Press Ctrl+C to stop the server" -ForegroundColor Yellow
Write-Host ""

python -m uvicorn app:app --reload --host 0.0.0.0 --port 8000

