

Write-Host "[INFO] Starting MolForge Backend Server..." -ForegroundColor Green
Write-Host ""

# Check if we're in the repo root
if (-not (Test-Path "backend\app.py")) {
    Write-Host "[ERROR] backend\app.py not found. Please run this script from the repo root" -ForegroundColor Red
    exit 1
}

# Activate virtual environment (check both backend/.venv and .venv)
if (Test-Path "backend\.venv\Scripts\Activate.ps1") {
    . backend\.venv\Scripts\Activate.ps1
    Write-Host "[OK] Virtual environment activated (backend/.venv)" -ForegroundColor Green
} elseif (Test-Path ".venv\Scripts\Activate.ps1") {
    . .venv\Scripts\Activate.ps1
    Write-Host "[OK] Virtual environment activated (.venv)" -ForegroundColor Green
} else {
    Write-Host "[WARN] No virtual environment found. Using system Python..." -ForegroundColor Yellow
    Write-Host "   Tip: Create one with: cd backend; python -m venv .venv" -ForegroundColor Yellow
}

# Ensure uvicorn is installed
python -c "import uvicorn" 2>$null
if ($LASTEXITCODE -ne 0) {
    Write-Host "[WARN] uvicorn not installed. Installing..." -ForegroundColor Yellow
    pip install uvicorn
}

# Ensure backend\.env exists
if (-not (Test-Path "backend\.env")) {
    if (Test-Path "backend\env.example") {
        Write-Host "[INFO] Creating backend/.env from env.example..." -ForegroundColor Yellow
        Copy-Item "backend\env.example" "backend\.env"
    } else {
        Write-Host "[WARN] backend\env.example missing; skipping .env creation" -ForegroundColor Yellow
    }
}

Write-Host ""
Write-Host "[INFO] Starting server on http://localhost:8000" -ForegroundColor Cyan
Write-Host "[INFO] API docs available at http://localhost:8000/docs" -ForegroundColor Cyan
Write-Host "[INFO] Health check: http://localhost:8000/health" -ForegroundColor Cyan
Write-Host ""
Write-Host "Press Ctrl+C to stop the server" -ForegroundColor Yellow
Write-Host ""

# Set PYTHONPATH to include repo root for backend imports
$env:PYTHONPATH = (Get-Location).Path

# Use --reload-dir to only watch backend directory and avoid node_modules issues
python -m uvicorn backend.app:app --reload --reload-dir backend --host 0.0.0.0 --port 8000



'''
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
.venv\Scripts\Activate.ps1

uvicorn backend.app:app --reload
'''
