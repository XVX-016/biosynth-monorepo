import { useEffect } from 'react';
import { useNavigate } from 'react-router-dom';
import { useAuthStore } from '../store/authStore';

/**
 * Hook that redirects to /login if user is not authenticated.
 * Useful for pages that need auth but aren't wrapped in ProtectedRoute.
 */
export function useRequireAuth() {
  const { user, loading, initialize, initialized } = useAuthStore();
  const navigate = useNavigate();

  useEffect(() => {
    if (!initialized) {
      initialize();
    }
  }, [initialized, initialize]);

  useEffect(() => {
    if (!loading && initialized && !user) {
      navigate('/login', { replace: true });
    }
  }, [user, loading, initialized, navigate]);

  return { user, loading: loading || !initialized };
}

