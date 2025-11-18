# Railway Build Fix - Monorepo Configuration

## Problem

Railway was trying to build all services (engine, frontend, backend) because it detected multiple Dockerfiles in the monorepo.

## Solution

I've configured Railway to:
1. ✅ Only build the backend service
2. ✅ Use the correct Dockerfile path (`backend/Dockerfile`)
3. ✅ Handle monorepo structure correctly

## Files Updated

### 1. `railway.json` (Root)
- Points to `backend/Dockerfile`
- Sets correct start command

### 2. `railway.toml` (Root)
- Alternative configuration format
- Same settings as railway.json

### 3. `backend/Dockerfile`
- Updated to work from monorepo root
- Copies files correctly: `COPY backend/requirements.txt`
- Sets PYTHONPATH for imports

### 4. `.railwayignore` (Root)
- Excludes frontend, engine, and other unnecessary files
- Reduces build size and time

## Railway Configuration Steps

### Option 1: Use railway.json (Recommended)

1. In Railway dashboard → Your Service → Settings
2. Go to **"Source"** section
3. Ensure **Root Directory** is set to: `.` (root of repo)
4. Railway will automatically use `railway.json`

### Option 2: Manual Configuration

1. In Railway dashboard → Your Service → Settings
2. Set **Root Directory**: `.` (root)
3. Set **Dockerfile Path**: `backend/Dockerfile`
4. Set **Start Command**: `cd backend && uvicorn app:app --host 0.0.0.0 --port $PORT`

## Verify Configuration

After updating, Railway should:
- ✅ Only build backend service
- ✅ Use backend/Dockerfile
- ✅ Start the FastAPI server correctly

## If Build Still Fails

### Check Railway Logs

1. Go to Railway dashboard → Your Service → Deployments
2. Click on the failed deployment
3. Check the build logs for errors

### Common Issues

**Issue: "Module not found: backend"**
- **Fix**: PYTHONPATH is set in Dockerfile, should work
- **Alternative**: Update imports to use relative paths

**Issue: "Cannot find requirements.txt"**
- **Fix**: Dockerfile now copies from `backend/requirements.txt`
- **Verify**: Check that requirements.txt exists in backend/

**Issue: "Port already in use"**
- **Fix**: Railway provides $PORT automatically
- **Verify**: Start command uses $PORT

**Issue: "CORS errors"**
- **Fix**: Set CORS_ORIGINS environment variable in Railway
- **Format**: `https://your-app.web.app,http://localhost:5173`

## Testing Locally

To test the Dockerfile locally:

```bash
# From repo root
docker build -f backend/Dockerfile -t biosynth-backend .
docker run -p 8000:8000 -e PORT=8000 biosynth-backend
```

## Next Steps

1. ✅ Push updated files to GitHub
2. ✅ Railway will auto-redeploy
3. ✅ Check build logs
4. ✅ Test deployed backend URL

