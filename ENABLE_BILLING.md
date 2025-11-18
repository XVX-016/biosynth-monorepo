# Enable Billing for Google Cloud

## ⚠️ Required Step

Cloud Run, Cloud Build, and Container Registry require billing to be enabled.

**Good News:** Google Cloud provides **$300 free credit** that lasts 90 days, and Cloud Run has a generous free tier!

## ✅ Enable Billing

### Step 1: Go to Billing

1. Open: https://console.cloud.google.com/billing?project=molforge-478621
2. Or: Cloud Console → Billing

### Step 2: Link Billing Account

1. Click **"Link a billing account"**
2. If you don't have one:
   - Click **"Create billing account"**
   - Enter payment method (credit card)
   - **You won't be charged** unless you exceed free tier

### Step 3: Link to Project

1. Select your billing account
2. Click **"Set account"**
3. Project `molforge-478621` is now linked

## 💰 Free Tier Limits

**Cloud Run Free Tier:**
- ✅ 2 million requests/month
- ✅ 360,000 GB-seconds
- ✅ 180,000 vCPU-seconds

**Cloud Build Free Tier:**
- ✅ 120 build-minutes/day

**You likely won't exceed these limits!**

## ✅ After Billing is Enabled

Come back and run:

```powershell
# Enable APIs
gcloud services enable cloudbuild.googleapis.com run.googleapis.com containerregistry.googleapis.com

# Then deploy
cd backend
.\deploy-cloud-run.ps1
```

---

## 🔒 Payment Protection

- **Free tier**: $300 credit (90 days)
- **Alerts**: Set up billing alerts to monitor usage
- **Budget**: Set a budget limit (e.g., $10/month)

**You can set a budget of $0 to prevent any charges!**

---

**Once billing is enabled, we can proceed with deployment! 🚀**

