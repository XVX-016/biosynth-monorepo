# Backend Deployment Guide

## Current Setup

The backend is currently configured for:
- ✅ **Local Development** (localhost:8000)
- ✅ **Docker Deployment** (Dockerfile included)
- ❌ **NOT configured for Firebase** (Firebase Hosting is for static frontend only)

## Important: Firebase Hosting Limitation

**Firebase Hosting only hosts static files** (HTML, CSS, JS). It cannot run Python/FastAPI backends.

For backend APIs, you need one of these options:

---

## Option 1: Firebase Cloud Functions (Serverless)

Convert your FastAPI app to Firebase Cloud Functions:

```bash
# Install Firebase CLI
npm install -g firebase-tools

# Initialize Firebase Functions
firebase init functions
```

**Pros:**
- Integrated with Firebase ecosystem
- Serverless (pay per use)
- Auto-scaling

**Cons:**
- Requires rewriting to Cloud Functions format
- Cold start latency
- Limited execution time (9 minutes max)

---

## Option 2: Google Cloud Run (Recommended for FastAPI)

Deploy your Docker container to Cloud Run:

```bash
# Build and push to Google Container Registry
gcloud builds submit --tag gcr.io/YOUR_PROJECT/biosynth-backend

# Deploy to Cloud Run
gcloud run deploy biosynth-backend \
  --image gcr.io/YOUR_PROJECT/biosynth-backend \
  --platform managed \
  --region us-central1 \
  --allow-unauthenticated
```

**Pros:**
- Works with existing Dockerfile
- Serverless (pay per request)
- Auto-scaling
- Can use Firebase Auth easily

**Cons:**
- Requires Google Cloud account
- Cold start latency

---

## Option 3: Railway.app (Easiest)

1. Push code to GitHub
2. Connect Railway to your repo
3. Railway auto-detects Dockerfile
4. Deploy!

**Pros:**
- Very easy setup
- Free tier available
- Auto-deploy from Git

**Cons:**
- Third-party service
- Limited free tier

---

## Option 4: Render.com

Similar to Railway, very easy deployment:

1. Connect GitHub repo
2. Select "Web Service"
3. Point to `backend/Dockerfile`
4. Set start command: `uvicorn app:app --host 0.0.0.0 --port $PORT`

**Pros:**
- Free tier
- Easy setup
- Auto-deploy

---

## Option 5: Keep Local + Use ngrok (Development Only)

For testing Firebase frontend with local backend:

```bash
# Install ngrok
# Run backend locally
uvicorn app:app --host 0.0.0.0 --port 8000

# In another terminal
ngrok http 8000
# Use the ngrok URL in frontend
```

**Pros:**
- Quick testing
- No deployment needed

**Cons:**
- Not for production
- Free tier has limitations

---

## Recommended Setup for Production

### Frontend: Firebase Hosting ✅
- Already configured in `frontend/firebase.json`
- Deploy: `firebase deploy --only hosting`

### Backend: Google Cloud Run ✅
- Uses your existing Dockerfile
- Integrates well with Firebase Auth
- Serverless and scalable

### Database: 
- **Development**: SQLite (current)
- **Production**: PostgreSQL (Cloud SQL) or Firebase Firestore

---

## Quick Start: Deploy to Cloud Run

1. **Install Google Cloud SDK**
   ```bash
   # Download from: https://cloud.google.com/sdk/docs/install
   ```

2. **Authenticate**
   ```bash
   gcloud auth login
   gcloud config set project YOUR_PROJECT_ID
   ```

3. **Build and Deploy**
   ```bash
   cd backend
   gcloud builds submit --tag gcr.io/YOUR_PROJECT/biosynth-backend
   gcloud run deploy biosynth-backend \
     --image gcr.io/YOUR_PROJECT/biosynth-backend \
     --platform managed \
     --region us-central1 \
     --allow-unauthenticated \
     --set-env-vars "DATABASE_URL=your-db-url"
   ```

4. **Update Frontend API URL**
   ```typescript
   // frontend/src/lib/api.ts
   const API_URL = import.meta.env.VITE_API_URL || 
     'https://biosynth-backend-xxxxx.run.app'
   ```

5. **Update CORS in Backend**
   ```python
   # backend/config.py
   CORS_ORIGINS: List[str] = [
       "https://your-firebase-app.web.app",
       "https://your-firebase-app.firebaseapp.com",
   ]
   ```

---

## Current Status

- ✅ Backend code is ready for deployment
- ✅ Dockerfile is configured
- ✅ Environment variables are set up
- ❌ Not connected to Firebase (by design - Firebase Hosting ≠ Backend hosting)
- ⚠️ Need to choose deployment platform (Cloud Run recommended)

