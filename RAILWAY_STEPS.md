# 🚀 Railway Deployment - Exact Steps

## Option 1: Root Directory = `backend` (RECOMMENDED - Simpler)

### Step 1: Set Root Directory
1. Railway Dashboard → `@biosynth/backend` → **Settings**
2. Click **"Add Root Directory"**
3. Enter: `backend`
4. Click **Save**

### Step 2: Update Dockerfile
1. Rename `backend/Dockerfile` to `backend/Dockerfile.old`
2. Rename `backend/Dockerfile.simple` to `backend/Dockerfile`
3. Commit and push:
   ```bash
   git add backend/Dockerfile backend/Dockerfile.simple
   git commit -m "Use simplified Dockerfile for Railway"
   git push
   ```

### Step 3: Delete Other Services
- Delete `@biosynth/frontend` service
- Delete `@biosynth/engine` service

### Step 4: Set Start Command
In Railway → Backend → Settings → Deploy:
- **Start Command**: `uvicorn app:app --host 0.0.0.0 --port $PORT`

### Step 5: Set Environment Variables
In Railway → Backend → Variables:
```
CORS_ORIGINS=https://your-app.web.app,http://localhost:5173
ENVIRONMENT=production
```

### Step 6: Redeploy
- Railway will auto-redeploy on push
- Or manually click **"Redeploy"**

---

## Option 2: Root Directory = `.` (Repo Root)

### Step 1: Set Root Directory
1. Railway Dashboard → `@biosynth/backend` → **Settings**
2. **Root Directory**: Leave empty or set to `.`
3. Click **Save**

### Step 2: Keep Current Dockerfile
- The current `backend/Dockerfile` works with repo root
- No changes needed

### Step 3: Set Dockerfile Path
In Railway → Backend → Settings → Build:
- **Dockerfile Path**: `backend/Dockerfile`

### Step 4: Set Start Command
In Railway → Backend → Settings → Deploy:
- **Start Command**: `cd backend && uvicorn app:app --host 0.0.0.0 --port $PORT`

### Step 5-6: Same as Option 1

---

## 🎯 Which Option to Choose?

**Option 1 (Root = `backend`)** is **RECOMMENDED** because:
- ✅ Simpler Dockerfile
- ✅ Cleaner build context
- ✅ Less confusion
- ✅ Standard Railway pattern

**Option 2 (Root = `.`)** if:
- You want to keep current Dockerfile
- You might add other services later

---

## 🔍 Check Logs After Deploy

1. Railway → Backend → **Logs** tab
2. Look for:
   - ✅ "Application startup complete" = Success!
   - ❌ Any error messages = Need to fix

---

## ✅ Success Indicators

- Backend service shows **green checkmark**
- Status: **"Active"** or **"Deployed"**
- Logs show: "Uvicorn running on http://0.0.0.0:XXXX"
- Health endpoint works: `curl https://your-url/health`

---

## 🆘 If Still Crashing

1. **Check Logs** - Share the error message
2. **Verify Root Directory** is set correctly
3. **Verify Start Command** uses `$PORT`
4. **Check Environment Variables** are set

