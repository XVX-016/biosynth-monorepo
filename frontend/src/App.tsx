import React from 'react';
import { Routes, Route } from 'react-router-dom';
import AppShell from './layouts/AppShell';
import Dashboard from './pages/Dashboard';
import Library from './pages/Library';
import Lab from './pages/Lab';
import Profile from './pages/Profile';
import Models from './pages/Models';
import Docs from './pages/Docs';
import AdminItems from './pages/admin/Items';

export default function App() {
	return (
		<AppShell>
			<Routes>
				<Route path="/" element={<Dashboard />} />
				<Route path="/lab" element={<Lab />} />
				<Route path="/library" element={<Library />} />
				<Route path="/models" element={<Models />} />
				<Route path="/docs" element={<Docs />} />
				<Route path="/profile" element={<Profile />} />
				<Route path="/admin/items" element={<AdminItems />} />
			</Routes>
		</AppShell>
	);
}
