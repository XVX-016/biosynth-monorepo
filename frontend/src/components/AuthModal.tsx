import React, { useState, FormEvent, useEffect } from 'react';
import { motion, AnimatePresence, useReducedMotion } from 'framer-motion';
import { useAuthStore } from '../store/authStore';
import Button from './ui/Button';
import MoleculeDecoration from './auth/MoleculeDecoration';
import {
  backdropVariants,
  modalVariants,
  tabContentVariants,
  errorShakeVariants,
} from './auth/authMotion';

type AuthModalTab = 'signin' | 'signup' | 'confirm-email';

export default function AuthModal() {
  const {
    authModalOpen,
    authModalTab,
    pendingEmail,
    signIn,
    signUp,
    resendConfirmationEmail,
    closeAuthModal,
    setAuthModalTab,
  } = useAuthStore();

  const [email, setEmail] = useState('');
  const [password, setPassword] = useState('');
  const [confirmPassword, setConfirmPassword] = useState('');
  const [error, setError] = useState<string | null>(null);
  const [isLoading, setIsLoading] = useState(false);
  const [resendLoading, setResendLoading] = useState(false);
  const [resendSuccess, setResendSuccess] = useState(false);
  const prefersReducedMotion = useReducedMotion();

  // Reset form when modal opens/closes or tab changes
  useEffect(() => {
    if (authModalOpen) {
      setError(null);
      setIsLoading(false);
      if (pendingEmail && authModalTab === 'confirm-email') {
        setEmail(pendingEmail);
      }
    } else {
      // Reset form when modal closes
      setEmail('');
      setPassword('');
      setConfirmPassword('');
      setError(null);
      setResendSuccess(false);
    }
  }, [authModalOpen, authModalTab, pendingEmail]);

  // Close on escape key
  useEffect(() => {
    if (!authModalOpen) return;

    const handleEsc = (e: KeyboardEvent) => {
      if (e.key === 'Escape') closeAuthModal();
    };
    window.addEventListener('keydown', handleEsc);
    return () => window.removeEventListener('keydown', handleEsc);
  }, [authModalOpen, closeAuthModal]);

  const handleSignIn = async (e: FormEvent) => {
    e.preventDefault();
    setError(null);
    
    if (isLoading) return;
    setIsLoading(true);

    try {
      const { error } = await signIn(email, password);
      if (error) {
        // Handle email not confirmed error
        if (error.message?.includes('email_not_confirmed') || error.message?.includes('Email not confirmed')) {
          setError('Please confirm your email before signing in.');
          setAuthModalTab('confirm-email');
        } else {
          setError(error.message || 'Failed to sign in');
        }
      }
      // Modal will close automatically via auth state change handler
    } catch (err: any) {
      setError(err.message || 'An unexpected error occurred');
    } finally {
      setIsLoading(false);
    }
  };

  const handleSignUp = async (e: FormEvent) => {
    e.preventDefault();
    setError(null);
    
    if (isLoading) return;

    // Validate passwords match
    if (password !== confirmPassword) {
      setError('Passwords do not match');
      return;
    }

    // Validate password length
    if (password.length < 6) {
      setError('Password must be at least 6 characters');
      return;
    }

    setIsLoading(true);

    try {
      const { error, requiresConfirmation } = await signUp(email, password);
      if (error) {
        setError(error.message || 'Failed to create account');
      } else if (requiresConfirmation) {
        // Switch to confirmation view (handled in store)
        setError(null);
      }
      // If no confirmation required, modal will close automatically via auth state change
    } catch (err: any) {
      setError(err.message || 'An unexpected error occurred');
    } finally {
      setIsLoading(false);
    }
  };

  const handleResendEmail = async () => {
    if (!email) return;
    
    setResendLoading(true);
    setResendSuccess(false);
    setError(null);

    try {
      const { error } = await resendConfirmationEmail(email);
      if (error) {
        setError(error.message || 'Failed to resend confirmation email');
      } else {
        setResendSuccess(true);
      }
    } catch (err: any) {
      setError(err.message || 'Failed to resend confirmation email');
    } finally {
      setResendLoading(false);
    }
  };

  if (!authModalOpen) return null;

  // Simplified variants for reduced motion
  const safeModalVariants = prefersReducedMotion
    ? { hidden: { opacity: 0 }, visible: { opacity: 1 }, exit: { opacity: 0 } }
    : modalVariants;
  const safeTabVariants = prefersReducedMotion
    ? { hidden: { opacity: 0 }, visible: { opacity: 1 }, exit: { opacity: 0 } }
    : tabContentVariants;

  return (
    <AnimatePresence>
      <motion.div
        variants={prefersReducedMotion ? {} : backdropVariants}
        initial="hidden"
        animate="visible"
        exit="exit"
        className="fixed inset-0 z-50 flex items-center justify-center bg-black/40 backdrop-blur-sm p-4"
        onClick={(e) => {
          if (e.target === e.currentTarget) closeAuthModal();
        }}
      >
        <motion.div
          variants={safeModalVariants}
          initial="hidden"
          animate="visible"
          exit="exit"
          className="bg-white rounded-2xl shadow-xl w-full max-w-md overflow-hidden relative"
          onClick={(e) => e.stopPropagation()}
        >
          {/* Header with tabs and close button */}
          <div className="relative px-8 pt-8 pb-6 border-b border-lightGrey flex items-center justify-between">
            {/* Tabs - hidden during email confirmation */}
            {authModalTab !== 'confirm-email' && (
              <div className="flex gap-6">
                <button
                  onClick={() => {
                    setAuthModalTab('signin');
                    setError(null);
                  }}
                  className={`px-4 py-2 font-medium transition-colors relative ${
                    authModalTab === 'signin'
                      ? 'text-black'
                      : 'text-darkGrey hover:text-black'
                  }`}
                >
                  Sign In
                  {authModalTab === 'signin' && (
                    <motion.div
                      layoutId="activeTab"
                      className="absolute bottom-0 left-0 right-0 h-0.5 bg-black -mb-6"
                      transition={{ duration: 0.2 }}
                    />
                  )}
                </button>
                <button
                  onClick={() => {
                    setAuthModalTab('signup');
                    setError(null);
                  }}
                  className={`px-4 py-2 font-medium transition-colors relative ${
                    authModalTab === 'signup'
                      ? 'text-black'
                      : 'text-darkGrey hover:text-black'
                  }`}
                >
                  Create Account
                  {authModalTab === 'signup' && (
                    <motion.div
                      layoutId="activeTab"
                      className="absolute bottom-0 left-0 right-0 h-0.5 bg-black -mb-6"
                      transition={{ duration: 0.2 }}
                    />
                  )}
                </button>
              </div>
            )}
            
            {/* Spacer when tabs are hidden (confirm-email view) */}
            {authModalTab === 'confirm-email' && <div />}

            {/* Close button */}
            <button
              onClick={closeAuthModal}
              className="p-2 text-darkGrey hover:text-black hover:bg-offwhite rounded-full transition-colors"
              aria-label="Close modal"
            >
              <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                <line x1="18" y1="6" x2="6" y2="18" />
                <line x1="6" y1="6" x2="18" y2="18" />
              </svg>
            </button>
          </div>

          <div className="relative overflow-hidden">
            {/* Background decoration */}
            <div className="absolute inset-0 flex items-center justify-center opacity-[0.03] pointer-events-none">
              <div className="w-full max-w-md">
                <MoleculeDecoration />
              </div>
            </div>

            {/* FORM */}
            <div className="relative p-8 max-w-md mx-auto">
              <AnimatePresence mode="wait">
                {/* Email Confirmation View */}
                {authModalTab === 'confirm-email' && (
                  <motion.div
                    key="confirm-email"
                    variants={safeTabVariants}
                    initial="hidden"
                    animate="visible"
                    exit="exit"
                    className="text-center"
                  >
                    <div className="text-5xl mb-4">✉️</div>
                    <h2 className="text-2xl font-bold text-black mb-2">Confirm your email</h2>
                    <p className="text-darkGrey mb-2">
                      We've sent a confirmation link to
                    </p>
                    <p className="text-black font-medium mb-6">{email}</p>
                    <p className="text-sm text-darkGrey mb-6">
                      Click the link in that email to activate your account.
                    </p>

                    {error && (
                      <motion.div
                        variants={prefersReducedMotion ? {} : errorShakeVariants}
                        animate="shake"
                        className="bg-red-50 border border-red-200 text-red-800 px-4 py-3 rounded-lg mb-4"
                      >
                        {error}
                      </motion.div>
                    )}

                    {resendSuccess && (
                      <motion.div
                        initial={{ opacity: 0, y: -10 }}
                        animate={{ opacity: 1, y: 0 }}
                        className="bg-green-50 border border-green-200 text-green-800 px-4 py-3 rounded-lg mb-4"
                      >
                        Confirmation email sent! Check your inbox.
                      </motion.div>
                    )}

                    <div className="space-y-3">
                      <Button
                        variant="primary"
                        onClick={handleResendEmail}
                        disabled={resendLoading}
                        className="w-full"
                      >
                        {resendLoading ? 'Sending...' : 'Resend email'}
                      </Button>
                      <Button
                        variant="secondary"
                        onClick={() => {
                          setAuthModalTab('signin');
                          setError(null);
                        }}
                        className="w-full"
                      >
                        Back to Sign In
                      </Button>
                    </div>

                    <p className="text-xs text-midGrey mt-4">
                      Didn't receive it? Check spam or wait a minute before resending.
                    </p>
                  </motion.div>
                )}

                {/* Sign In Form */}
                {authModalTab === 'signin' && (
                  <motion.div
                    key="signin"
                    variants={safeTabVariants}
                    initial="hidden"
                    animate="visible"
                    exit="exit"
                  >
                    <h2 className="text-2xl font-bold text-black mb-2">Sign in to MolForge</h2>
                    <form onSubmit={handleSignIn} className="mt-6 space-y-4">
                      {error && (
                        <motion.div
                          variants={prefersReducedMotion ? {} : errorShakeVariants}
                          animate="shake"
                          className="bg-red-50 border border-red-200 text-red-800 px-4 py-3 rounded-lg"
                        >
                          {error}
                        </motion.div>
                      )}
                      <div>
                        <label htmlFor="signin-email" className="block text-sm font-medium text-darkGrey mb-2">
                          Email address
                        </label>
                        <input
                          id="signin-email"
                          type="email"
                          autoComplete="email"
                          required
                          value={email}
                          onChange={(e) => setEmail(e.target.value)}
                          disabled={isLoading}
                          className="appearance-none relative block w-full px-4 py-3 border border-lightGrey placeholder-darkGrey text-black rounded-lg focus:outline-none focus:ring-2 focus:ring-black focus:border-transparent disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
                          placeholder="you@example.com"
                        />
                      </div>
                      <div>
                        <label htmlFor="signin-password" className="block text-sm font-medium text-darkGrey mb-2">
                          Password
                        </label>
                        <input
                          id="signin-password"
                          type="password"
                          autoComplete="current-password"
                          required
                          value={password}
                          onChange={(e) => setPassword(e.target.value)}
                          disabled={isLoading}
                          className="appearance-none relative block w-full px-4 py-3 border border-lightGrey placeholder-darkGrey text-black rounded-lg focus:outline-none focus:ring-2 focus:ring-black focus:border-transparent disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
                          placeholder="••••••••"
                        />
                      </div>
                      <Button
                        type="submit"
                        variant="primary"
                        className="w-full mt-4 h-11"
                        disabled={isLoading}
                      >
                        {isLoading ? 'Signing in...' : 'Sign in'}
                      </Button>
                    </form>
                  </motion.div>
                )}

                {/* Sign Up Form */}
                {authModalTab === 'signup' && (
                  <motion.div
                    key="signup"
                    variants={safeTabVariants}
                    initial="hidden"
                    animate="visible"
                    exit="exit"
                  >
                    <h2 className="text-2xl font-bold text-black mb-2">Create your MolForge account</h2>
                    <p className="text-sm text-darkGrey mb-6">
                      Your account lets you save molecules, favorites, and lab sessions.
                    </p>
                    <form onSubmit={handleSignUp} className="space-y-4">
                      {error && (
                        <motion.div
                          variants={prefersReducedMotion ? {} : errorShakeVariants}
                          animate="shake"
                          className="bg-red-50 border border-red-200 text-red-800 px-4 py-3 rounded-lg"
                        >
                          {error}
                        </motion.div>
                      )}
                      <div>
                        <label htmlFor="signup-email" className="block text-sm font-medium text-darkGrey mb-2">
                          Email address
                        </label>
                        <input
                          id="signup-email"
                          type="email"
                          autoComplete="email"
                          required
                          value={email}
                          onChange={(e) => setEmail(e.target.value)}
                          disabled={isLoading}
                          className="appearance-none relative block w-full px-4 py-3 border border-lightGrey placeholder-darkGrey text-black rounded-lg focus:outline-none focus:ring-2 focus:ring-black focus:border-transparent disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
                          placeholder="you@example.com"
                        />
                      </div>
                      <div>
                        <label htmlFor="signup-password" className="block text-sm font-medium text-darkGrey mb-2">
                          Password
                        </label>
                        <input
                          id="signup-password"
                          type="password"
                          autoComplete="new-password"
                          required
                          value={password}
                          onChange={(e) => setPassword(e.target.value)}
                          disabled={isLoading}
                          className="appearance-none relative block w-full px-4 py-3 border border-lightGrey placeholder-darkGrey text-black rounded-lg focus:outline-none focus:ring-2 focus:ring-black focus:border-transparent disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
                          placeholder="••••••••"
                        />
                        <p className="mt-1 text-xs text-darkGrey">Must be at least 6 characters</p>
                      </div>
                      <div>
                        <label htmlFor="signup-confirm" className="block text-sm font-medium text-darkGrey mb-2">
                          Confirm Password
                        </label>
                        <input
                          id="signup-confirm"
                          type="password"
                          autoComplete="new-password"
                          required
                          value={confirmPassword}
                          onChange={(e) => setConfirmPassword(e.target.value)}
                          disabled={isLoading}
                          className="appearance-none relative block w-full px-4 py-3 border border-lightGrey placeholder-darkGrey text-black rounded-lg focus:outline-none focus:ring-2 focus:ring-black focus:border-transparent disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
                          placeholder="••••••••"
                        />
                      </div>
                      <Button
                        type="submit"
                        variant="primary"
                        className="w-full mt-4 h-11"
                        disabled={isLoading}
                      >
                        {isLoading ? 'Creating account...' : 'Create account'}
                      </Button>
                    </form>
                    <p className="text-xs text-midGrey mt-4 text-center">
                      You'll need to confirm your email before signing in.
                    </p>
                  </motion.div>
                )}
              </AnimatePresence>
            </div>
          </div>
        </motion.div>
      </motion.div>
    </AnimatePresence>
  );
}
