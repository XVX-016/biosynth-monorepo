# Next Steps - Cloud Run Deployment

## ✅ Completed

- [x] Google Cloud SDK installed
- [x] Project created: `molforge-478621`
- [x] Project set in gcloud
- [x] Authenticated: `tanmmay2005@gmail.com`

## ⏳ Required: Enable Billing

**You need to enable billing before proceeding.**

### Quick Steps:

1. **Go to:** https://console.cloud.google.com/billing?project=molforge-478621

2. **Link billing account:**
   - Click "Link a billing account"
   - Create one if needed (requires credit card)
   - **You won't be charged** - free tier includes $300 credit

3. **Set budget (optional but recommended):**
   - Set budget to $0 or $10/month
   - Get alerts if usage exceeds

## 🚀 After Billing is Enabled

### Step 1: Enable APIs

```powershell
# Make sure gcloud is in PATH (add this to your session)
$env:PATH += ";C:\Users\tanmm\AppData\Local\Google\Cloud SDK\google-cloud-sdk\bin"

# Enable required APIs
gcloud services enable cloudbuild.googleapis.com run.googleapis.com containerregistry.googleapis.com
```

### Step 2: Deploy Backend

**Option A: Using PowerShell Script**

```powershell
cd backend
.\deploy-cloud-run.ps1
```

**Option B: Manual Deploy**

```powershell
# Build image
gcloud builds submit --tag gcr.io/molforge-478621/biosynth-backend ./backend

# Deploy to Cloud Run
gcloud run deploy biosynth-backend `
  --image gcr.io/molforge-478621/biosynth-backend `
  --platform managed `
  --region us-central1 `
  --allow-unauthenticated `
  --port 8080 `
  --memory 2Gi `
  --cpu 2 `
  --timeout 300 `
  --max-instances 10 `
  --set-env-vars "PORT=8080,ENVIRONMENT=production"
```

### Step 3: Get Service URL

```powershell
gcloud run services describe biosynth-backend --region us-central1 --format 'value(status.url)'
```

### Step 4: Set CORS

```powershell
gcloud run services update biosynth-backend `
  --region us-central1 `
  --update-env-vars "CORS_ORIGINS=https://your-firebase-app.web.app,https://your-firebase-app.firebaseapp.com,http://localhost:5173"
```

### Step 5: Test

```powershell
$SERVICE_URL = gcloud run services describe biosynth-backend --region us-central1 --format 'value(status.url)'
curl $SERVICE_URL/health
# Should return: {"status":"healthy"}
```

---

## 📝 Quick Reference

**Project Info:**
- Project ID: `molforge-478621`
- Project Number: `986946918424`
- Account: `tanmmay2005@gmail.com`

**gcloud Path:**
```powershell
C:\Users\tanmm\AppData\Local\Google\Cloud SDK\google-cloud-sdk\bin\gcloud.cmd
```

**Add to PATH (for this session):**
```powershell
$env:PATH += ";C:\Users\tanmm\AppData\Local\Google\Cloud SDK\google-cloud-sdk\bin"
```

---

## 💡 Tip: Add gcloud to PATH Permanently

1. Press `Win + X` → System → Advanced system settings
2. Environment Variables → User variables → Path → Edit
3. Add: `C:\Users\tanmm\AppData\Local\Google\Cloud SDK\google-cloud-sdk\bin`
4. Restart PowerShell

---

**Enable billing, then we can deploy! 🚀**

