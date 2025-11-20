/**
 * SupabaseStatus Component
 * 
 * A small component that shows Supabase connection status.
 * Can be added to the navbar or dashboard.
 */

import React, { useState, useEffect } from 'react';
import { supabase, isSupabaseConfigured } from '../supabase';

export default function SupabaseStatus() {
  const [isConfigured, setIsConfigured] = useState(false);
  const [isAuthenticated, setIsAuthenticated] = useState(false);
  const [userId, setUserId] = useState<string | null>(null);

  useEffect(() => {
    setIsConfigured(isSupabaseConfigured());
    
    if (supabase) {
      // Check current session
      supabase.auth.getSession().then(({ data: { session } }) => {
        setIsAuthenticated(!!session);
        setUserId(session?.user?.id || null);
      });

      // Listen for auth changes
      const {
        data: { subscription },
      } = supabase.auth.onAuthStateChange((_event, session) => {
        setIsAuthenticated(!!session);
        setUserId(session?.user?.id || null);
      });

      return () => subscription.unsubscribe();
    }
  }, []);

  if (!isConfigured) {
    return (
      <div className="flex items-center gap-2 text-sm">
        <span className="w-2 h-2 rounded-full bg-red-500"></span>
        <span className="text-midGrey">Supabase: Not Configured</span>
      </div>
    );
  }

  return (
    <div className="flex items-center gap-2 text-sm">
      <span className={`w-2 h-2 rounded-full ${isAuthenticated ? 'bg-green-500' : 'bg-yellow-500'}`}></span>
      <span className="text-darkGrey">
        Supabase: {isAuthenticated ? 'Connected' : 'Not Authenticated'}
      </span>
      {userId && (
        <span className="text-xs text-midGrey">({userId.substring(0, 8)}...)</span>
      )}
    </div>
  );
}

