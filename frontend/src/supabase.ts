import { createClient, type SupabaseClient } from '@supabase/supabase-js';

// Supabase configuration from environment variables
const supabaseUrl = import.meta.env.VITE_SUPABASE_URL;
const supabaseAnonKey = import.meta.env.VITE_SUPABASE_ANON_KEY;

// Debug: Log what Vite is actually reading (only in dev mode)
if (import.meta.env.DEV) {
  console.log('🔍 Debug - Environment variables from Vite:');
  console.log('  VITE_SUPABASE_URL:', supabaseUrl ? `${supabaseUrl.substring(0, 30)}...` : 'undefined');
  console.log('  VITE_SUPABASE_ANON_KEY:', supabaseAnonKey ? `${supabaseAnonKey.substring(0, 30)}...` : 'undefined');
  console.log('  All VITE_ vars:', Object.keys(import.meta.env).filter(k => k.startsWith('VITE_')));
}

// Validate Supabase configuration
function validateSupabaseConfig() {
  if (!supabaseUrl || !supabaseAnonKey) {
    console.warn(
      '⚠️ Supabase configuration incomplete. Missing:',
      !supabaseUrl ? 'VITE_SUPABASE_URL' : '',
      !supabaseAnonKey ? 'VITE_SUPABASE_ANON_KEY' : '',
      '\nPlease create a .env file with your Supabase credentials.'
    );
    return false;
  }

  // Check for placeholder values
  const hasPlaceholders = supabaseUrl.includes('your-') || 
                          supabaseAnonKey.includes('your-') ||
                          supabaseUrl.includes('placeholder') ||
                          supabaseAnonKey.includes('placeholder');
  
  if (hasPlaceholders) {
    console.warn('⚠️ Supabase configuration contains placeholder values.');
    console.warn('📝 Please update your .env file with actual Supabase credentials:');
    console.warn('   1. Go to https://supabase.com/dashboard');
    console.warn('   2. Select your project → Settings → API');
    console.warn('   3. Copy Project URL and anon/public key');
    console.warn('   4. Update frontend/.env with real values');
    console.warn('   5. Restart dev server');
    return false;
  }

  // Validate URL format
  if (!supabaseUrl.startsWith('https://') || !supabaseUrl.includes('.supabase.co')) {
    console.warn('⚠️ VITE_SUPABASE_URL appears invalid. Should be: https://xxxxx.supabase.co');
    return false;
  }

  // Validate key format (anon keys start with eyJ, NOT sb_secret_)
  if (supabaseAnonKey.startsWith('sb_secret_')) {
    console.error('❌ ERROR: You are using a SECRET KEY instead of ANON KEY!');
    console.error('   Secret keys (sb_secret_...) should NEVER be used in client-side code!');
    console.error('   They are only for server-side use and will expose your database!');
    console.error('');
    console.error('📝 How to fix:');
    console.error('   1. Go to Supabase Dashboard → Settings → API');
    console.error('   2. Find "anon public" key (NOT "service_role" or "secret")');
    console.error('   3. The anon key starts with "eyJ" (it\'s a JWT token)');
    console.error('   4. Update VITE_SUPABASE_ANON_KEY in .env with the anon key');
    console.error('   5. Restart dev server');
    return false;
  }
  
  if (!supabaseAnonKey.startsWith('eyJ')) {
    console.warn('⚠️ VITE_SUPABASE_ANON_KEY appears invalid. Should start with "eyJ"');
    console.warn('   Make sure you\'re using the "anon public" key, not the "service_role" or "secret" key');
    return false;
  }

  return true;
}

let supabase: SupabaseClient | null = null;

try {
  const isValid = validateSupabaseConfig();
  
  if (isValid && supabaseUrl && supabaseAnonKey) {
    supabase = createClient(supabaseUrl, supabaseAnonKey);
    console.log('✅ Supabase initialized successfully');
  } else {
    console.warn('⚠️ Supabase not initialized - configuration missing');
  }
} catch (error) {
  console.error('❌ Failed to initialize Supabase:', error);
}

// Export with null check
export { supabase };
export const isSupabaseConfigured = () => supabase !== null;

