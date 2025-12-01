# Stop any process using port 8000
$port = 8000
$processes = Get-NetTCPConnection -LocalPort $port -ErrorAction SilentlyContinue | Select-Object -ExpandProperty OwningProcess -Unique

if ($processes) {
    Write-Host "Stopping processes on port $port..."
    foreach ($pid in $processes) {
        try {
            Stop-Process -Id $pid -Force
            Write-Host "  Stopped process $pid"
        } catch {
            Write-Host "  Could not stop process $pid: $_"
        }
    }
} else {
    Write-Host "No processes found on port $port"
}

