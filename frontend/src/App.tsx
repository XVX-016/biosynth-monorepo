import { Routes, Route, useNavigate, useLocation } from 'react-router-dom';
import { useEffect } from 'react';
import AppShell from './layouts/AppShell';
import Dashboard from './pages/Dashboard';
import LibraryPage from './pages/LibraryPage';
import LabV2Page from './components/LabV2/LabV2Page';
import Profile from './pages/Profile';
import Models from './pages/Models';
import Docs from './pages/Docs';
import AdminItems from './pages/admin/Items';
import SupabaseTest from './pages/SupabaseTest';
import SeedLibrary from './pages/SeedLibrary';
import PublicLibrary from './pages/PublicLibrary';
import Phase10Dashboard from './pages/Phase10Dashboard';
import ProtectedRoute from './components/ProtectedRoute';
import AuthModal from './components/AuthModal';
import { useAuthStore } from './store/authStore';

// Component to handle /login and /signup routes
function AuthRouteHandler() {
	const location = useLocation();
	const navigate = useNavigate();
	const { openAuthModal } = useAuthStore();

	useEffect(() => {
		if (location.pathname === '/login') {
			openAuthModal('signin');
			navigate('/', { replace: true });
		} else if (location.pathname === '/signup') {
			openAuthModal('signup');
			navigate('/', { replace: true });
		}
	}, [location.pathname, navigate, openAuthModal]);

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
							<Route path="/models" element={<Models />} />
							<Route path="/docs" element={<Docs />} />
							<Route path="/phase10" element={<Phase10Dashboard />} />
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
