# 🚨 Railway Fix - Step by Step (URGENT)

## Current Status
- ✅ Backend: Built but **CRASHED** (needs root directory)
- ❌ Frontend: Failed (we don't need it on Railway)
- ❌ Engine: Failed (we don't need it on Railway)

## 🔧 IMMEDIATE FIXES

### Step 1: Set Root Directory (CRITICAL)

1. In Railway Dashboard → Click on **`@biosynth/backend`** service
2. Go to **Settings** tab (you're already there)
3. Find **"Add Root Directory"** section
4. Click **"Add Root Directory"**
5. Enter: `backend`
6. Click **Save** or **Update**

**This tells Railway to build from the `backend/` folder, not the repo root.**

### Step 2: Delete Unnecessary Services

Since we only need backend on Railway:

1. Go back to **Architecture** view
2. Click on **`@biosynth/frontend`** service
3. Go to **Settings** → Scroll down → Click **"Delete Service"**
4. Repeat for **`@biosynth/engine`** service

**Why?** Frontend should deploy to Firebase, not Railway. Engine is a package, not a service.

### Step 3: Check Why Backend Crashed

1. Click on **`@biosynth/backend`** service
2. Go to **Logs** tab
3. Look for error messages

**Common crash reasons:**
- ❌ "ModuleNotFoundError: No module named 'backend'"
- ❌ "Cannot find requirements.txt"
- ❌ Port binding issues
- ❌ Database connection errors

### Step 4: Update Start Command

1. In **Settings** → **Deploy** section
2. Set **Start Command** to:
   ```
   cd backend && uvicorn app:app --host 0.0.0.0 --port $PORT
   ```

### Step 5: Redeploy

1. Go to **Deployments** tab
2. Click **"Redeploy"** button
3. Or push a new commit to trigger auto-deploy

---

## 📋 Complete Configuration Checklist

### Backend Service Settings:

✅ **Source:**
- Root Directory: `backend`
- Branch: `master` (or your main branch)

✅ **Build:**
- Dockerfile Path: `Dockerfile` (relative to root directory)
- Build Command: (leave empty, uses Dockerfile)

✅ **Deploy:**
- Start Command: `cd backend && uvicorn app:app --host 0.0.0.0 --port $PORT`

✅ **Variables:**
- `CORS_ORIGINS`: `https://your-firebase-app.web.app,http://localhost:5173`
- `ENVIRONMENT`: `production`
- `DATABASE_URL`: `sqlite:///./biosynth.db` (or Railway PostgreSQL URL)

---

## 🔍 Debugging the Crash

### Check Logs

1. Railway Dashboard → Backend Service → **Logs** tab
2. Look for the most recent error
3. Common issues:

**Issue: "ModuleNotFoundError: No module named 'backend'"**
- **Fix**: Root Directory must be set to `backend`
- **Fix**: PYTHONPATH is set in Dockerfile, should work

**Issue: "Cannot find app:app"**
- **Fix**: Start command should be: `cd backend && uvicorn app:app ...`

**Issue: "Port already in use"**
- **Fix**: Use `$PORT` environment variable (already in start command)

**Issue: "Database connection failed"**
- **Fix**: Set `DATABASE_URL` in Variables
- **Fix**: If using SQLite, ensure path is correct

---

## 🎯 Quick Action Items

1. ✅ **Set Root Directory to `backend`** (MOST IMPORTANT)
2. ✅ **Delete frontend and engine services**
3. ✅ **Check logs for crash reason**
4. ✅ **Update start command if needed**
5. ✅ **Redeploy backend**

---

## 📝 After Fix

Once backend is running:

1. Copy the Railway public URL
2. Update frontend `.env` or `api.ts`:
   ```typescript
   const API_URL = 'https://your-railway-url.up.railway.app'
   ```
3. Deploy frontend to Firebase (not Railway)
4. Test the full stack

---

## 🆘 Still Crashing?

Share the error logs and I'll help debug!

