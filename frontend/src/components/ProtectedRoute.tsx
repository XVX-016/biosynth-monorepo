import React, { useEffect } from 'react';
import { Navigate, useLocation } from 'react-router-dom';
import { useAuthStore } from '../store/authStore';

interface ProtectedRouteProps {
  children: React.ReactElement;
}

/**
 * ProtectedRoute component that guards routes requiring authentication.
 * Redirects to /login if user is not authenticated.
 */
export default function ProtectedRoute({ children }: ProtectedRouteProps) {
  const { user, loading, initialize, initialized } = useAuthStore();
  const location = useLocation();

  // Initialize auth store on mount if not already initialized
  useEffect(() => {
    if (!initialized) {
      initialize();
    }
  }, [initialized, initialize]);

  // Show nothing while loading to prevent flash
  if (loading || !initialized) {
    return (
      <div className="min-h-screen flex items-center justify-center bg-offwhite">
        <div className="text-darkGrey">Loading...</div>
      </div>
    );
  }

  // Redirect to login if not authenticated, preserving the intended destination
  if (!user) {
    return <Navigate to="/login" state={{ from: location }} replace />;
  }

  return children;
}

