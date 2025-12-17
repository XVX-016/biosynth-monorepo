import { useEffect } from 'react';
import { useNavigate } from 'react-router-dom';
import { useAuthStore } from '../store/authStore';

/**
 * Component that handles navigation after successful authentication.
 * Watches for intendedDestination in auth store and navigates accordingly.
 */
export default function AuthNavigationHandler() {
  const navigate = useNavigate();
  const { user, intendedDestination } = useAuthStore();

  useEffect(() => {
    // When user is authenticated and there's an intended destination, navigate to it
    if (user && intendedDestination) {
      // Clear the intended destination and navigate
      const destination = intendedDestination;
      useAuthStore.setState({ intendedDestination: null });
      navigate(destination, { replace: true });
    }
  }, [user, intendedDestination, navigate]);

  return null;
}








