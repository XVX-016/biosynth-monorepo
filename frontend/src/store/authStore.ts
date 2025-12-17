import { create } from 'zustand';
import { devtools } from 'zustand/middleware';
import type { Session, User } from '@supabase/supabase-js';
import { supabase } from '../supabase';

type AuthModalTab = 'signin' | 'signup' | 'confirm-email';

interface AuthState {
  user: User | null;
  session: Session | null;
  loading: boolean;
  initialized: boolean;
  
  // Modal state
  authModalOpen: boolean;
  authModalTab: AuthModalTab;
  pendingAction: (() => Promise<void>) | null;
  pendingEmail: string | null;
  
  // Actions
  signIn: (email: string, password: string) => Promise<{ error: any }>;
  signUp: (email: string, password: string) => Promise<{ error: any; requiresConfirmation?: boolean }>;
  signOut: () => Promise<void>;
  initialize: () => Promise<void>;
  
  // Modal control
  openAuthModal: (tab?: AuthModalTab) => void;
  closeAuthModal: () => void;
  setAuthModalTab: (tab: AuthModalTab) => void;
  
  // Pending action handling
  setPendingAction: (action: (() => Promise<void>) | null) => void;
  runPendingAction: () => Promise<void>;
  
  // Email confirmation
  resendConfirmationEmail: (email: string) => Promise<{ error: any }>;
}

// Store the subscription outside the store to prevent multiple listeners
let authStateSubscription: { unsubscribe: () => void } | null = null;

export const useAuthStore = create<AuthState>()(
  devtools(
    (set, get) => ({
      user: null,
      session: null,
      loading: true,
      initialized: false,
      authModalOpen: false,
      authModalTab: 'signin',
      pendingAction: null,
      pendingEmail: null,

      openAuthModal: (tab: AuthModalTab = 'signin') => {
        set({ authModalOpen: true, authModalTab: tab });
      },

      closeAuthModal: () => {
        set({ authModalOpen: false, pendingEmail: null });
      },

      setAuthModalTab: (tab: AuthModalTab) => {
        set({ authModalTab: tab });
      },

      setPendingAction: (action: (() => Promise<void>) | null) => {
        set({ pendingAction: action });
      },

      runPendingAction: async () => {
        const action = get().pendingAction;
        if (action) {
          set({ pendingAction: null });
          try {
            await action();
          } catch (error) {
            console.error('Error running pending action:', error);
          }
        }
      },

      resendConfirmationEmail: async (email: string) => {
        if (!supabase) {
          return { error: { message: 'Supabase is not configured' } };
        }
        try {
          const { error } = await supabase.auth.resend({
            type: 'signup',
            email,
          });
          return { error };
        } catch (error: any) {
          return { error };
        }
      },

      initialize: async () => {
        const state = get();
        if (state.initialized || !supabase) {
          return;
        }

        // Prevent multiple simultaneous initializations
        set({ loading: true });

        try {
          // Set up auth state change listener first (only once)
          if (!authStateSubscription) {
            const { data: { subscription } } = supabase.auth.onAuthStateChange((event, session) => {
              const newUser = session?.user ?? null;
              set({ 
                user: newUser, 
                session,
                loading: false 
              });
              
              // Handle successful sign in/sign up
              if ((event === 'SIGNED_IN' || event === 'TOKEN_REFRESHED') && newUser) {
                // Close modal and run pending action
                const { runPendingAction, closeAuthModal } = get();
                closeAuthModal();
                runPendingAction();
              }
            });
            authStateSubscription = subscription;
          }

          // Get initial session
          const { data: { session }, error } = await supabase.auth.getSession();
          
          if (error) {
            console.error('Error getting session:', error);
            set({ user: null, session: null, loading: false, initialized: true });
            return;
          }

          set({ 
            user: session?.user ?? null, 
            session, 
            loading: false,
            initialized: true
          });
        } catch (error) {
          console.error('Error initializing auth:', error);
          set({ user: null, session: null, loading: false, initialized: true });
        }
      },

      signIn: async (email: string, password: string) => {
        if (!supabase) {
          return { error: { message: 'Supabase is not configured' } };
        }

        try {
          const { error } = await supabase.auth.signInWithPassword({ email, password });
          // State will be updated via onAuthStateChange listener
          return { error };
        } catch (error: any) {
          return { error };
        }
      },

      signUp: async (email: string, password: string) => {
        if (!supabase) {
          return { error: { message: 'Supabase is not configured' } };
        }

        try {
          const { data, error } = await supabase.auth.signUp({ email, password });
          
          // Check if email confirmation is required
          const requiresConfirmation = !data.user?.email_confirmed_at && !error;
          
          if (requiresConfirmation) {
            set({ authModalTab: 'confirm-email', pendingEmail: email });
          }
          
          // State will be updated via onAuthStateChange listener if user is immediately authenticated
          return { error, requiresConfirmation };
        } catch (error: any) {
          return { error };
        }
      },

      signOut: async () => {
        if (supabase) {
          await supabase.auth.signOut();
          // State will be updated via onAuthStateChange listener
        }
      },
    }),
    { name: 'AuthStore' }
  )
);

