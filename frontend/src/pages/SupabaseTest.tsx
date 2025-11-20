/**
 * Supabase Test Page
 * 
 * This page helps verify that Supabase is properly configured and working.
 * Access it at /supabase-test
 */

import React, { useState, useEffect } from 'react';
import { supabase, isSupabaseConfigured } from '../supabase';
import { saveMolecule, listMolecules, deleteMolecule } from '../lib/supabaseMoleculeStore';

interface TestResult {
  name: string;
  status: 'pending' | 'success' | 'error';
  message: string;
}

export default function SupabaseTest() {
  const [results, setResults] = useState<TestResult[]>([]);
  const [isRunning, setIsRunning] = useState(false);
  const [userId, setUserId] = useState<string | null>(null);
  const [testMoleculeId, setTestMoleculeId] = useState<string | null>(null);

  useEffect(() => {
    if (!supabase) {
      return;
    }
    
    // Check current session
    supabase.auth.getSession().then(({ data: { session } }) => {
      setUserId(session?.user?.id || null);
    });

    // Listen for auth changes
    const {
      data: { subscription },
    } = supabase.auth.onAuthStateChange((_event, session) => {
      setUserId(session?.user?.id || null);
    });

    return () => subscription.unsubscribe();
  }, []);

  const addResult = (name: string, status: 'success' | 'error', message: string) => {
    setResults((prev) => [...prev, { name, status, message }]);
  };

  const clearResults = () => {
    setResults([]);
    setTestMoleculeId(null);
  };

  const testSupabaseConfig = () => {
    clearResults();
    addResult('Supabase Configuration', isSupabaseConfigured() ? 'success' : 'error', 
      isSupabaseConfigured() 
        ? 'Supabase is properly configured' 
        : 'Supabase configuration is missing or invalid. Check your .env file.');
  };

  const testDatabaseConnection = async () => {
    if (!isSupabaseConfigured() || !supabase) {
      addResult('Database Connection', 'error', 'Supabase not initialized. Please configure Supabase in your .env file.');
      return;
    }

    try {
      // Test connection by querying a system table
      const { data, error } = await supabase.from('molecules').select('count').limit(1);
      
      if (error && error.code !== 'PGRST116') {
        // PGRST116 means table doesn't exist yet, which is okay
        addResult('Database Connection', 'error', `Error: ${error.message}`);
      } else {
        addResult('Database Connection', 'success', 'Successfully connected to Supabase database');
      }
    } catch (error: any) {
      addResult('Database Connection', 'error', `Error: ${error.message}`);
    }
  };

  const testAuthentication = async () => {
    if (!isSupabaseConfigured() || !supabase) {
      addResult('Authentication', 'error', 'Supabase not initialized. Please configure Supabase in your .env file.');
      return;
    }

    try {
      if (!userId) {
        // Sign in anonymously
        const { data, error } = await supabase.auth.signInAnonymously();
        if (error) throw error;
        addResult('Anonymous Sign In', 'success', `Signed in as: ${data.user?.id}`);
        setUserId(data.user?.id || null);
      } else {
        addResult('Authentication', 'success', `Already signed in as: ${userId}`);
      }
    } catch (error: any) {
      addResult('Authentication', 'error', `Error: ${error.message}`);
    }
  };

  const testMoleculeStore = async () => {
    if (!isSupabaseConfigured()) {
      addResult('Molecule Store', 'error', 'Supabase not configured. Please update your .env file.');
      return;
    }
    if (!userId) {
      addResult('Molecule Store', 'error', 'Please sign in first');
      return;
    }

    try {
      // Test save
      const moleculeId = await saveMolecule(userId, {
        name: 'Test Molecule',
        smiles: 'CCO',
        formula: 'C2H6O',
      });
      setTestMoleculeId(moleculeId);
      addResult('Save Molecule', 'success', `Saved molecule with ID: ${moleculeId}`);

      // Test list
      const molecules = await listMolecules(userId);
      addResult('List Molecules', 'success', `Found ${molecules.length} molecule(s)`);

      // Test delete
      if (moleculeId) {
        await deleteMolecule(userId, moleculeId);
        addResult('Delete Molecule', 'success', 'Successfully deleted test molecule');
        setTestMoleculeId(null);
      }
    } catch (error: any) {
      addResult('Molecule Store', 'error', `Error: ${error.message}`);
    }
  };

  const runAllTests = async () => {
    setIsRunning(true);
    clearResults();

    // Test 1: Configuration
    testSupabaseConfig();
    await new Promise((resolve) => setTimeout(resolve, 500));

    // Test 2: Database connection
    await testDatabaseConnection();
    await new Promise((resolve) => setTimeout(resolve, 500));

    // Test 3: Authentication
    await testAuthentication();
    await new Promise((resolve) => setTimeout(resolve, 500));

    // Test 4: Molecule store (requires auth)
    if (userId || supabase) {
      await testMoleculeStore();
    }

    setIsRunning(false);
  };

  const handleSignOut = async () => {
    if (!supabase) {
      addResult('Sign Out', 'error', 'Supabase not initialized');
      return;
    }
    try {
      await supabase.auth.signOut();
      setUserId(null);
      addResult('Sign Out', 'success', 'Signed out successfully');
    } catch (error: any) {
      addResult('Sign Out', 'error', `Error: ${error.message}`);
    }
  };

  return (
    <div className="p-8 max-w-4xl mx-auto space-y-6">
      <div>
        <h1 className="text-3xl font-bold text-black mb-2">Supabase Test Page</h1>
        <p className="text-darkGrey">
          Use this page to verify that Supabase is properly configured and working.
        </p>
      </div>

      {/* Status */}
      <div className="bg-white rounded-xl shadow-neon border border-lightGrey p-4">
        <h2 className="text-xl font-semibold mb-3">Status</h2>
        <div className="space-y-2">
          <div className="flex items-center gap-2">
            <span className="font-medium">Supabase Configured:</span>
            <span className={isSupabaseConfigured() ? 'text-green-600' : 'text-red-600'}>
              {isSupabaseConfigured() ? '✅ Yes' : '❌ No'}
            </span>
          </div>
          <div className="flex items-center gap-2">
            <span className="font-medium">Authenticated:</span>
            <span className={userId ? 'text-green-600' : 'text-gray-600'}>
              {userId ? `✅ Yes (${userId.substring(0, 8)}...)` : '❌ No'}
            </span>
          </div>
          {userId && (
            <button
              onClick={handleSignOut}
              className="btn-secondary text-sm px-3 py-1 mt-2"
            >
              Sign Out
            </button>
          )}
        </div>
      </div>

      {/* Test Controls */}
      <div className="bg-white rounded-xl shadow-neon border border-lightGrey p-4">
        <h2 className="text-xl font-semibold mb-3">Tests</h2>
        <div className="flex flex-wrap gap-2">
          <button
            onClick={runAllTests}
            disabled={isRunning || !isSupabaseConfigured()}
            className="btn-primary px-4 py-2"
          >
            {isRunning ? 'Running Tests...' : 'Run All Tests'}
          </button>
          <button
            onClick={testSupabaseConfig}
            disabled={isRunning}
            className="btn-secondary px-4 py-2"
          >
            Test Config
          </button>
          <button
            onClick={testDatabaseConnection}
            disabled={isRunning || !isSupabaseConfigured()}
            className="btn-secondary px-4 py-2"
          >
            Test Database
          </button>
          <button
            onClick={testAuthentication}
            disabled={isRunning || !isSupabaseConfigured()}
            className="btn-secondary px-4 py-2"
          >
            Test Auth
          </button>
          <button
            onClick={testMoleculeStore}
            disabled={isRunning || !userId}
            className="btn-secondary px-4 py-2"
          >
            Test Molecule Store
          </button>
          <button
            onClick={clearResults}
            disabled={isRunning}
            className="btn-secondary px-4 py-2"
          >
            Clear Results
          </button>
        </div>
      </div>

      {/* Results */}
      {results.length > 0 && (
        <div className="bg-white rounded-xl shadow-neon border border-lightGrey p-4">
          <h2 className="text-xl font-semibold mb-3">Test Results</h2>
          <div className="space-y-2">
            {results.map((result, index) => (
              <div
                key={index}
                className={`p-3 rounded-lg border ${
                  result.status === 'success'
                    ? 'bg-green-50 border-green-200'
                    : result.status === 'error'
                    ? 'bg-red-50 border-red-200'
                    : 'bg-gray-50 border-gray-200'
                }`}
              >
                <div className="flex items-center gap-2 mb-1">
                  <span className="font-medium">{result.name}</span>
                  <span>
                    {result.status === 'success' && '✅'}
                    {result.status === 'error' && '❌'}
                    {result.status === 'pending' && '⏳'}
                  </span>
                </div>
                <p className="text-sm text-darkGrey">{result.message}</p>
              </div>
            ))}
          </div>
        </div>
      )}

      {/* Instructions */}
      <div className="bg-offwhite rounded-xl border border-lightGrey p-4">
        <h2 className="text-xl font-semibold mb-3">Setup Instructions</h2>
        {!isSupabaseConfigured() && (
          <div className="mb-4 p-3 bg-yellow-50 border border-yellow-200 rounded-lg">
            <p className="font-semibold text-yellow-800 mb-2">⚠️ Supabase Not Configured</p>
            <p className="text-sm text-yellow-700">
              Your <code className="bg-yellow-100 px-1 rounded">.env</code> file contains placeholder values.
              Please replace them with your actual Supabase credentials.
            </p>
          </div>
        )}
        <ol className="list-decimal list-inside space-y-2 text-sm text-darkGrey">
          <li>Create a Supabase project at <a href="https://supabase.com" target="_blank" rel="noopener noreferrer" className="text-blue-600 underline">Supabase</a></li>
          <li>Go to Project Settings → API</li>
          <li>Copy your Project URL and anon/public key</li>
          <li>Update <code className="bg-lightGrey px-1 rounded">frontend/.env</code> with:
            <ul className="list-disc list-inside ml-4 mt-1 space-y-1">
              <li><code>VITE_SUPABASE_URL</code> - Your project URL</li>
              <li><code>VITE_SUPABASE_ANON_KEY</code> - Your anon/public key</li>
            </ul>
          </li>
          <li>Create the molecules table in Supabase SQL Editor (see SQL below)</li>
          <li>Enable Row Level Security (RLS) policies</li>
          <li>Restart your dev server after updating <code className="bg-lightGrey px-1 rounded">.env</code></li>
          <li>Run the tests above to verify everything works</li>
        </ol>
        <div className="mt-4 p-3 bg-gray-50 rounded-lg">
          <p className="font-semibold mb-2">SQL to create molecules table:</p>
          <pre className="text-xs overflow-x-auto bg-gray-100 p-2 rounded">
{`CREATE TABLE molecules (
  id UUID DEFAULT gen_random_uuid() PRIMARY KEY,
  name TEXT NOT NULL,
  smiles TEXT,
  formula TEXT,
  json_graph TEXT,
  properties TEXT,
  thumbnail_b64 TEXT,
  user_id UUID NOT NULL REFERENCES auth.users(id),
  created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
);

-- Enable RLS
ALTER TABLE molecules ENABLE ROW LEVEL SECURITY;

-- Policy: Users can only see their own molecules
CREATE POLICY "Users can view own molecules"
  ON molecules FOR SELECT
  USING (auth.uid() = user_id);

-- Policy: Users can insert their own molecules
CREATE POLICY "Users can insert own molecules"
  ON molecules FOR INSERT
  WITH CHECK (auth.uid() = user_id);

-- Policy: Users can update their own molecules
CREATE POLICY "Users can update own molecules"
  ON molecules FOR UPDATE
  USING (auth.uid() = user_id);

-- Policy: Users can delete their own molecules
CREATE POLICY "Users can delete own molecules"
  ON molecules FOR DELETE
  USING (auth.uid() = user_id);`}
          </pre>
        </div>
      </div>
    </div>
  );
}

