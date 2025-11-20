/**
 * Quick script to verify Supabase configuration
 * Run with: node scripts/verify-supabase.js
 */

import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const envPath = path.join(__dirname, '..', '.env');

console.log('🔍 Verifying Supabase Configuration...\n');

// Check if .env exists
if (!fs.existsSync(envPath)) {
  console.error('❌ .env file not found at:', envPath);
  console.log('💡 Create it from .env.example');
  process.exit(1);
}

// Read .env file
const envContent = fs.readFileSync(envPath, 'utf8');

// Check for required variables
const requiredVars = {
  'VITE_SUPABASE_URL': envContent.match(/VITE_SUPABASE_URL=(.+)/)?.[1]?.trim(),
  'VITE_SUPABASE_ANON_KEY': envContent.match(/VITE_SUPABASE_ANON_KEY=(.+)/)?.[1]?.trim(),
};

console.log('📋 Environment Variables:');
let allValid = true;

Object.entries(requiredVars).forEach(([key, value]) => {
  if (!value) {
    console.log(`  ❌ ${key}: Missing`);
    allValid = false;
  } else if (value.includes('your-') || value.includes('placeholder')) {
    console.log(`  ⚠️  ${key}: Contains placeholder value`);
    console.log(`     Current: ${value.substring(0, 30)}...`);
    allValid = false;
  } else {
    const displayValue = key.includes('KEY') 
      ? `${value.substring(0, 20)}...` 
      : value;
    console.log(`  ✅ ${key}: ${displayValue}`);
  }
});

console.log('\n📦 Package Check:');
const packageJsonPath = path.join(__dirname, '..', 'package.json');
if (fs.existsSync(packageJsonPath)) {
  const packageJson = JSON.parse(fs.readFileSync(packageJsonPath, 'utf8'));
  if (packageJson.dependencies && packageJson.dependencies['@supabase/supabase-js']) {
    console.log(`  ✅ @supabase/supabase-js: ${packageJson.dependencies['@supabase/supabase-js']}`);
  } else {
    console.log('  ❌ @supabase/supabase-js: Not found in dependencies');
    console.log('  💡 Run: npm install');
    allValid = false;
  }
}

console.log('\n📁 File Check:');
const requiredFiles = [
  'src/supabase.ts',
  'src/lib/supabaseMoleculeStore.ts',
  'src/components/SupabaseStatus.tsx',
  'src/pages/SupabaseTest.tsx',
];

requiredFiles.forEach(file => {
  const filePath = path.join(__dirname, '..', file);
  if (fs.existsSync(filePath)) {
    console.log(`  ✅ ${file}`);
  } else {
    console.log(`  ❌ ${file}: Missing`);
    allValid = false;
  }
});

console.log('\n' + '='.repeat(50));
if (allValid) {
  console.log('✅ All checks passed! Supabase is configured correctly.');
  console.log('\n💡 Next steps:');
  console.log('   1. Make sure you\'ve created the molecules table in Supabase');
  console.log('   2. Run: npm run dev');
  console.log('   3. Visit: http://localhost:5173/supabase-test');
} else {
  console.log('⚠️  Some issues found. Please fix them before proceeding.');
}
console.log('='.repeat(50));

process.exit(allValid ? 0 : 1);

