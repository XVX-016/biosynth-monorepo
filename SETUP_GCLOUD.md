# Google Cloud SDK Setup for Windows

## Issue: `gcloud` command not found

The Google Cloud SDK is installed but not in your PATH.

## ✅ Quick Fix

### Option 1: Add to PATH (Recommended)

1. **Find gcloud installation:**
   - Usually at: `C:\Users\YOUR_USERNAME\AppData\Local\Google\Cloud SDK\google-cloud-sdk\bin`
   - Or: `C:\Program Files\Google\Cloud SDK\google-cloud-sdk\bin`

2. **Add to PATH:**
   - Press `Win + X` → System → Advanced system settings
   - Click "Environment Variables"
   - Under "User variables", select "Path" → Edit
   - Click "New" → Add the path above
   - Click OK on all windows

3. **Restart PowerShell/Terminal**

### Option 2: Use Full Path

Instead of `gcloud`, use the full path:

```powershell
# Find where gcloud is installed
Get-ChildItem -Path "$env:LOCALAPPDATA\Google" -Recurse -Filter "gcloud.cmd" -ErrorAction SilentlyContinue
Get-ChildItem -Path "$env:ProgramFiles\Google" -Recurse -Filter "gcloud.cmd" -ErrorAction SilentlyContinue
```

Then use the full path, or create an alias:

```powershell
# In current session only
$env:PATH += ";C:\Users\YOUR_USERNAME\AppData\Local\Google\Cloud SDK\google-cloud-sdk\bin"
```

### Option 3: Reinstall with PATH Option

1. Download Google Cloud SDK installer
2. During installation, check "Add to PATH"
3. Restart terminal

---

## 🚀 After gcloud is Working

### Step 1: Set Project

```powershell
gcloud config set project molforge-478621
```

### Step 2: Authenticate

```powershell
gcloud auth login
gcloud auth application-default login
```

### Step 3: Enable APIs

```powershell
gcloud services enable cloudbuild.googleapis.com
gcloud services enable run.googleapis.com
gcloud services enable containerregistry.googleapis.com
```

### Step 4: Deploy

```powershell
# From repo root
cd backend
.\deploy-cloud-run.ps1
```

---

## 🔍 Verify Installation

After adding to PATH, restart PowerShell and test:

```powershell
gcloud --version
# Should show: Google Cloud SDK version
```

---

## 📝 Alternative: Use Cloud Console

If CLI setup is difficult, you can deploy via Cloud Console:

1. Go to: https://console.cloud.google.com/run
2. Select project: `molforge-478621`
3. Click "Create Service"
4. Upload Dockerfile or use Cloud Build
5. Configure and deploy

---

**Once gcloud is in PATH, you can proceed with deployment!**

