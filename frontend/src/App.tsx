import { Routes, Route, useNavigate, useLocation } from 'react-router-dom';
import { useEffect } from 'react';
import AppShell from './layouts/AppShell';
import Dashboard from './pages/Dashboard';
import LibraryPage from './pages/LibraryPage';
import LabV2Page from './components/LabV2/LabV2Page';
import Profile from './pages/Profile';
import Docs from './pages/Docs';
import AdminItems from './pages/admin/Items';
import SupabaseTest from './pages/SupabaseTest';
import SeedLibrary from './pages/SeedLibrary';
import PublicLibrary from './pages/PublicLibrary';
import StudioPage from './pages/StudioPage';
import ProtectedRoute from './components/ProtectedRoute';
import AuthModal from './components/AuthModal';
import AuthNavigationHandler from './components/AuthNavigationHandler';
import { useAuthStore } from './store/authStore';

// Component to handle /login and /signup routes
function AuthRouteHandler() {
	const location = useLocation();
	const navigate = useNavigate();
	const { openAuthModal } = useAuthStore();

	useEffect(() => {
		if (location.pathname === '/login') {
			// Extract intended destination from state, fallback to null
			const intendedDestination = (location.state as { from?: { pathname: string } })?.from?.pathname || null;
			openAuthModal('signin', intendedDestination);
			// Navigate to home but preserve the intended destination in store
			navigate('/', { replace: true });
		} else if (location.pathname === '/signup') {
			// Extract intended destination from state, fallback to null
			const intendedDestination = (location.state as { from?: { pathname: string } })?.from?.pathname || null;
			openAuthModal('signup', intendedDestination);
			// Navigate to home but preserve the intended destination in store
			navigate('/', { replace: true });
		}
	}, [location.pathname, location.state, navigate, openAuthModal]);

	return null;
}

// Main Application Router
export default function App() {
	return (
		<>
			{/* Auth Modal - accessible from anywhere */}
			<AuthModal />
			{/* Handle /login and /signup routes */}
			<AuthRouteHandler />
			{/* Handle post-login navigation to intended destination */}
			<AuthNavigationHandler />

			<Routes>
				{/* Lab route - fullscreen, no AppShell, PUBLIC */}
				<Route path="/lab" element={<LabV2Page />} />

				{/* All other routes - wrapped in AppShell */}
				<Route path="/*" element={
					<AppShell>
						<Routes>
							<Route path="/" element={<Dashboard />} />
							<Route path="/library" element={<LibraryPage />} />
							<Route path="/library/public" element={<PublicLibrary />} />
							<Route path="/docs" element={<Docs />} />
							<Route path="/studio" element={<StudioPage />} />
							{/* Authenticated-only routes */}
							<Route
								path="/profile"
								element={
									<ProtectedRoute>
										<Profile />
									</ProtectedRoute>
								}
							/>
							<Route
								path="/my-molecules"
								element={
									<ProtectedRoute>
										<LibraryPage />
									</ProtectedRoute>
								}
							/>
							<Route path="/admin/items" element={<AdminItems />} />
							<Route path="/supabase-test" element={<SupabaseTest />} />
							<Route path="/seed-library" element={<SeedLibrary />} />
						</Routes>
					</AppShell>
				} />
			</Routes>
		</>
	);
}
