# Deploy to Cloud Run - Step by Step

## ✅ Project Setup Complete

- **Project ID**: `molforge-478621`
- **Project Number**: `986946918424`

## 🚀 Next Steps

### Step 1: Verify Project is Set

```bash
gcloud config set project molforge-478621
gcloud config get-value project
# Should show: molforge-478621
```

### Step 2: Enable Required APIs

```bash
gcloud services enable cloudbuild.googleapis.com
gcloud services enable run.googleapis.com
gcloud services enable containerregistry.googleapis.com
```

### Step 3: Authenticate (if not already)

```bash
gcloud auth login
gcloud auth application-default login
```

### Step 4: Deploy Backend

**Option A: Using PowerShell Script (Easiest)**

```powershell
cd backend
.\deploy-cloud-run.ps1
```

**Option B: Manual Deploy**

```bash
# From repo root
gcloud builds submit --tag gcr.io/molforge-478621/biosynth-backend ./backend

gcloud run deploy biosynth-backend \
  --image gcr.io/molforge-478621/biosynth-backend \
  --platform managed \
  --region us-central1 \
  --allow-unauthenticated \
  --port 8080 \
  --memory 2Gi \
  --cpu 2 \
  --timeout 300 \
  --max-instances 10 \
  --set-env-vars "PORT=8080,ENVIRONMENT=production"
```

### Step 5: Get Your Service URL

```bash
gcloud run services describe biosynth-backend \
  --region us-central1 \
  --format 'value(status.url)'
```

### Step 6: Set CORS Origins

```bash
gcloud run services update biosynth-backend \
  --region us-central1 \
  --update-env-vars "CORS_ORIGINS=https://your-firebase-app.web.app,https://your-firebase-app.firebaseapp.com,http://localhost:5173"
```

### Step 7: Test Deployment

```bash
# Get URL
SERVICE_URL=$(gcloud run services describe biosynth-backend --region us-central1 --format 'value(status.url)')

# Test
curl $SERVICE_URL/health
# Should return: {"status":"healthy"}
```

---

## 📝 Quick Command Reference

```bash
# Set project
gcloud config set project molforge-478621

# Enable APIs
gcloud services enable cloudbuild.googleapis.com run.googleapis.com containerregistry.googleapis.com

# Build and deploy
gcloud builds submit --tag gcr.io/molforge-478621/biosynth-backend ./backend
gcloud run deploy biosynth-backend --image gcr.io/molforge-478621/biosynth-backend --platform managed --region us-central1 --allow-unauthenticated --port 8080 --memory 2Gi --cpu 2

# Get URL
gcloud run services describe biosynth-backend --region us-central1 --format 'value(status.url)'
```

---

## ✅ Checklist

- [ ] Project set: `molforge-478621`
- [ ] APIs enabled
- [ ] Authenticated
- [ ] Backend deployed
- [ ] Service URL obtained
- [ ] CORS_ORIGINS set
- [ ] Health endpoint tested
- [ ] Frontend API URL updated

---

**Ready to deploy! Run the commands above. 🚀**

