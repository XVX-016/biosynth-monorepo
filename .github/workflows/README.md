# GitHub Actions CI/CD Pipeline

This directory contains the CI/CD workflows for the BioSynth AI monorepo.

## Workflows

### 1. `ci.yml` - Main CI Pipeline

Runs on every push to `main`/`master` and on pull requests. Performs:

- ✅ TypeScript type checking (Frontend + Engine)
- ✅ ESLint linting (Frontend)
- ✅ Vitest unit tests (Frontend + Engine)
- ✅ Build verification (Frontend + Engine)
- ✅ Artifact uploads

**Features:**
- Uses npm workspaces for monorepo support
- Caches node_modules and Vite build cache
- Fast fail-fast on errors
- Uploads build artifacts for deployment

### 2. `deploy.yml` - Automated Deployment

Deploys to Vercel when code is pushed to `main`/`master`.

**Setup Required:**
1. Add these secrets to your GitHub repository:
   - `VERCEL_TOKEN` - Your Vercel API token
   - `VERCEL_ORG_ID` - Your Vercel organization ID
   - `VERCEL_PROJECT_ID` - Your Vercel project ID

2. Get these from:
   - Vercel Dashboard → Settings → Tokens
   - Vercel Dashboard → Project Settings → General

**Alternative Deployments:**
- For Netlify: Use `netlify-cli` action
- For AWS S3: Use `aws-cli` action
- For Cloudflare Pages: Use `cloudflare/pages-action`

### 3. `pr-comments.yml` - PR Status Bot

Automatically comments on pull requests with:
- Test results summary
- Type checking status
- Linting status
- Pass/fail counts

Runs on PR open, update, and reopen events.

## Setup Instructions

### 1. Ensure ESLint is Installed

If you haven't already, install ESLint dependencies in `frontend/`:

```bash
cd frontend
npm install -D eslint @eslint/js globals eslint-plugin-react-hooks eslint-plugin-react-refresh typescript-eslint
```

### 2. Update README Badges

Edit `README.md` and replace `YOUR-ORG` with your actual GitHub organization/username:

```markdown
![CI](https://github.com/YOUR-ORG/biosynth-monorepo/actions/workflows/ci.yml/badge.svg)
```

### 3. Configure Vercel Secrets (for deployment)

1. Go to your GitHub repository
2. Settings → Secrets and variables → Actions
3. Add the three Vercel secrets mentioned above

### 4. Test the Pipeline

1. Push to a branch
2. Create a pull request
3. Check the Actions tab to see the pipeline run

## Customization

### Adding More Checks

To add additional checks (e.g., backend tests), add steps to `ci.yml`:

```yaml
- name: Run Backend Tests
  run: cd backend && pytest
```

### Changing Node Version

Update the `node-version` in all workflow files:

```yaml
node-version: 20  # Change to 18, 22, etc.
```

### Adjusting Cache Keys

Modify cache keys in `ci.yml` to invalidate when needed:

```yaml
key: vite-${{ runner.os }}-${{ hashFiles('**/package-lock.json') }}-${{ hashFiles('**/*.ts') }}
```

## Troubleshooting

### Tests Fail in CI but Pass Locally

- Check Node.js version matches (CI uses Node 20)
- Ensure all dependencies are in `package.json` (not just `package-lock.json`)
- Verify test environment matches (jsdom, etc.)

### Linting Fails

- Ensure ESLint dependencies are installed
- Check `eslint.config.js` syntax
- Verify file patterns match your project structure

### Build Fails

- Check TypeScript errors first (`tsc --noEmit`)
- Verify all imports resolve correctly
- Check for missing dependencies

### Deployment Fails

- Verify Vercel secrets are set correctly
- Check Vercel project settings
- Ensure build output directory matches Vercel config

## Performance Tips

1. **Use caching** - Already configured for node_modules and Vite
2. **Parallel jobs** - Split frontend/backend into separate jobs if needed
3. **Conditional runs** - Only run tests for changed packages (advanced)

## Security

- Never commit secrets to the repository
- Use GitHub Secrets for all sensitive data
- Review third-party actions before using them
- Keep actions updated to latest versions




