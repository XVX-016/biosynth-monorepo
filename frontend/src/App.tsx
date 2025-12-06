
import { Routes, Route } from 'react-router-dom';
import AppShell from './layouts/AppShell';
import Dashboard from './pages/Dashboard';
import LibraryPage from './pages/LibraryPage';
import LabPage from "./pages/lab/Lab";
import Profile from './pages/Profile';
import Models from './pages/Models';
import Docs from './pages/Docs';
import AdminItems from './pages/admin/Items';
import SupabaseTest from './pages/SupabaseTest';
import SeedLibrary from './pages/SeedLibrary';
import PublicLibrary from './pages/PublicLibrary';
import Phase10Dashboard from './pages/Phase10Dashboard';

export default function App() {
	return (
		<Routes>
			{/* Lab route - fullscreen, no AppShell */}
			<Route path="/lab" element={<LabPage />} />

			{/* All other routes - wrapped in AppShell */}
			<Route path="/*" element={
				<AppShell>
					<Routes>
						<Route path="/" element={<Dashboard />} />
						<Route path="/library" element={<LibraryPage />} />
						<Route path="/library/public" element={<PublicLibrary />} />
						<Route path="/models" element={<Models />} />
						<Route path="/docs" element={<Docs />} />
						<Route path="/profile" element={<Profile />} />
						<Route path="/admin/items" element={<AdminItems />} />
						<Route path="/supabase-test" element={<SupabaseTest />} />
						<Route path="/seed-library" element={<SeedLibrary />} />
						<Route path="/phase10" element={<Phase10Dashboard />} />
					</Routes>
				</AppShell>
			} />
		</Routes>
	);
}
