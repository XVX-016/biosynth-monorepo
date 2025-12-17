import { useAuthStore } from '../store/authStore';

/**
 * Helper function to require authentication before performing an action.
 * If user is not authenticated, opens auth modal and stores the action to retry after login.
 * 
 * @param action The action to perform (must be async)
 * @returns Promise that resolves when action completes or rejects if auth fails
 */
export async function requireAuth(action: () => Promise<void>): Promise<void> {
  const { user, openAuthModal, setPendingAction } = useAuthStore.getState();

  if (!user) {
    // Store the action to retry after successful login
    setPendingAction(action);
    // Open auth modal
    openAuthModal('signin');
    // Return a promise that will be resolved when auth completes
    // The actual action will run via runPendingAction in authStore
    return Promise.resolve();
  }

  // User is authenticated, run the action immediately
  return action();
}

