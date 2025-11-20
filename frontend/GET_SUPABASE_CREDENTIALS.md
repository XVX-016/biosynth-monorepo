# 🔑 How to Get Your Supabase Credentials

## Quick Steps

### 1. Go to Supabase Dashboard

1. Visit [Supabase Dashboard](https://supabase.com/dashboard)
2. Sign in or create an account
3. Select your project (or create a new one)

### 2. Get Your Credentials

1. In your Supabase project, click **Settings** (gear icon) in the left sidebar
2. Click **API** in the settings menu
3. You'll see two important values:

   **Project URL:**
   - Format: `https://xxxxx.supabase.co`
   - Copy this value

   **anon/public key:**
   - Starts with `eyJ...`
   - This is your public/anonymous key
   - Copy this value

### 3. Update Your .env File

Open `frontend/.env` and replace the placeholder values:

**BEFORE (Wrong):**
```env
VITE_SUPABASE_URL=https://your-project-id.supabase.co
VITE_SUPABASE_ANON_KEY=your-anon-key-here
```

**AFTER (Right):**
```env
VITE_SUPABASE_URL=https://abcdefghijklmnop.supabase.co
VITE_SUPABASE_ANON_KEY=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6ImFiY2RlZmdoaWprbG1ub3AiLCJyb2xlIjoiYW5vbiIsImlhdCI6MTYzODU2Nzg5MCwiZXhwIjoxOTU0MTQzODkwfQ.example
```

**Important:**
- ✅ Remove ALL quotes (no `"` around values)
- ✅ No spaces around the `=` sign
- ✅ Use exact values from Supabase Dashboard
- ✅ URL should start with `https://` and end with `.supabase.co`
- ✅ Key should start with `eyJ`

### 4. Restart Dev Server

**CRITICAL:** After updating `.env`:

1. Stop dev server: Press `Ctrl+C`
2. Start again: `npm run dev`
3. Check browser console for: `✅ Supabase initialized successfully`

## ✅ Verification

After updating and restarting:

1. **Browser Console** (F12):
   - Should see: `✅ Supabase initialized successfully`
   - No warnings about placeholder values

2. **Test Page**:
   - Visit: http://localhost:5173/supabase-test
   - Should show: `Supabase Configured: ✅ Yes`

## 🚨 Common Mistakes

### ❌ Wrong: Still has placeholders
```env
VITE_SUPABASE_URL=https://your-project-id.supabase.co  # WRONG!
```

### ✅ Right: Actual value
```env
VITE_SUPABASE_URL=https://abcdefghijklmnop.supabase.co  # RIGHT!
```

### ❌ Wrong: Has quotes
```env
VITE_SUPABASE_URL="https://..."  # WRONG!
```

### ✅ Right: No quotes
```env
VITE_SUPABASE_URL=https://...  # RIGHT!
```

### ❌ Wrong: Wrong format
```env
VITE_SUPABASE_URL=supabase.com/project  # WRONG!
```

### ✅ Right: Correct format
```env
VITE_SUPABASE_URL=https://xxxxx.supabase.co  # RIGHT!
```

## 📝 Example .env File

Here's what a correctly filled `.env` file looks like:

```env
# Supabase Configuration
VITE_SUPABASE_URL=https://abcdefghijklmnop.supabase.co
VITE_SUPABASE_ANON_KEY=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6ImFiY2RlZmdoaWprbG1ub3AiLCJyb2xlIjoiYW5vbiIsImlhdCI6MTYzODU2Nzg5MCwiZXhwIjoxOTU0MTQzODkwfQ.example
```

**Note:** These are example values. Use your actual values from Supabase Dashboard.

## 🔗 Quick Links

- [Supabase Dashboard](https://supabase.com/dashboard)
- [Supabase Docs](https://supabase.com/docs)
- [API Settings](https://supabase.com/dashboard/project/_/settings/api)

---

**Remember:** Replace ALL placeholder values with your actual Supabase credentials!

