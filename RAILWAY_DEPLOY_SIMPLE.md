# Railway Deployment - Simple Steps

## ✅ No Billing Required!

Railway doesn't require billing for the free tier.

## 🚀 Deploy in 3 Steps

### Step 1: Go to Railway

1. Visit: https://railway.app
2. Sign in with GitHub
3. Click **"New Project"**
4. Select **"Deploy from GitHub repo"**
5. Choose: `XVX-016/biosynth-monorepo`

### Step 2: Configure Service

1. Railway will create services automatically
2. **Delete** `frontend` and `engine` services (we only need backend)
3. Click on **`@biosynth/backend`** service
4. Go to **Settings** → **Deployment**
5. Set **Root Directory**: `backend`
6. **Save**

### Step 3: Set Environment Variables

1. Go to **Variables** tab
2. Add:
   ```
   CORS_ORIGINS=https://your-firebase-app.web.app,http://localhost:5173
   ENVIRONMENT=production
   ```

### Step 4: Deploy

Railway will automatically:
- Build your Docker image
- Deploy to Cloud Run
- Provide a public URL

### Step 5: Get Your URL

1. Railway Dashboard → Backend Service
2. Click on the service
3. Copy the **Public Domain** (e.g., `biosynth-backend-production.up.railway.app`)

### Step 6: Update Frontend

In `frontend/src/lib/api.ts`:

```typescript
const API_URL = import.meta.env.VITE_API_URL || 
  'https://your-railway-url.up.railway.app'
```

---

## ✅ That's It!

No billing, no credit card, just deploy!

---

## 🎯 Your Dockerfile is Already Optimized

- ✅ Uses RDKit base image (~700MB)
- ✅ CPU-only PyTorch (~120MB instead of 800MB)
- ✅ Proper caching layers
- ✅ Ready for Railway

---

**Deploy to Railway now - it's the easiest option! 🚀**

