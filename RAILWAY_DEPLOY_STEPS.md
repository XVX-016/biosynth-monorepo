# Railway Deployment - Step by Step Fix

## 🚨 Current Issue

Railway is trying to build all services (engine, frontend, backend) and failing.

## ✅ Solution Applied

I've updated the configuration to:
1. Only build backend service
2. Use correct Dockerfile path for monorepo
3. Handle Python imports correctly

## 📋 Steps to Fix in Railway Dashboard

### Step 1: Update Service Settings

1. Go to Railway dashboard
2. Click on your **backend service**
3. Go to **Settings** → **Source**

### Step 2: Configure Root Directory

Set **Root Directory** to: `.` (root of repository)

**OR** if Railway has a "Service Root" option, leave it empty/default.

### Step 3: Configure Dockerfile Path

In **Settings** → **Build**:
- **Dockerfile Path**: `backend/Dockerfile`
- **Build Command**: (leave empty, uses Dockerfile)
- **Start Command**: `cd backend && uvicorn app:app --host 0.0.0.0 --port $PORT`

### Step 4: Redeploy

1. Go to **Deployments** tab
2. Click **"Redeploy"** or push new commit to trigger rebuild
3. Railway will use the updated configuration

## 🔍 Verify Configuration

After updating, check:

1. **Build Logs** should show:
   ```
   Building from Dockerfile: backend/Dockerfile
   ```

2. **No errors** about:
   - "Cannot find requirements.txt" ✅
   - "Module not found: backend" ✅
   - Multiple services building ❌

3. **Successful build** with:
   ```
   Successfully built ...
   ```

## 📝 Files Changed

### Root Level (New)
- ✅ `railway.json` - Railway configuration
- ✅ `railway.toml` - Alternative config format
- ✅ `.railwayignore` - Exclude unnecessary files

### Backend
- ✅ `backend/Dockerfile` - Updated for monorepo structure

## 🚀 Quick Fix Commands

If you have Railway CLI:

```bash
# From repo root
railway link
railway service
railway variables set CORS_ORIGINS="https://your-app.web.app,http://localhost:5173"
railway up
```

## ⚙️ Environment Variables

Make sure these are set in Railway:

```env
CORS_ORIGINS=https://your-firebase-app.web.app,https://your-firebase-app.firebaseapp.com,http://localhost:5173
ENVIRONMENT=production
DATABASE_URL=sqlite:///./biosynth.db
```

## 🧪 Test After Deployment

```bash
# Test health endpoint
curl https://your-railway-url.up.railway.app/health

# Should return: {"status":"healthy"}
```

## ❌ If Still Failing

### Check Build Logs

1. Railway Dashboard → Service → Deployments
2. Click on latest deployment
3. Check **Build Logs** tab

### Common Errors & Fixes

**Error: "COPY failed: file not found"**
- ✅ Fixed: Dockerfile now uses `COPY backend/requirements.txt`

**Error: "ModuleNotFoundError: No module named 'backend'"**
- ✅ Fixed: PYTHONPATH=/app is set in Dockerfile

**Error: "Multiple services detected"**
- ✅ Fixed: railway.json specifies only backend

**Error: "Port already in use"**
- ✅ Fixed: Using $PORT environment variable

## 📞 Next Steps

1. ✅ Update Railway settings (see above)
2. ✅ Push changes to GitHub (if not already)
3. ✅ Trigger redeploy in Railway
4. ✅ Check build logs
5. ✅ Test deployed URL
6. ✅ Update frontend API URL

## 💡 Pro Tip

If Railway still tries to build multiple services:
1. Delete other services in Railway dashboard
2. Keep only the backend service
3. Or create separate Railway projects for each service

