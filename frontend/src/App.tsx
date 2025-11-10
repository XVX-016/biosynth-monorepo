import React from 'react';
import { Routes, Route } from 'react-router-dom';
import AppShell from './layouts/AppShell';
import Dashboard from './pages/Dashboard';
import Library from './pages/Library';
import Lab from './pages/Lab';
import Profile from './pages/Profile';

export default function App() {
	return (
		<AppShell>
			<Routes>
				<Route path="/" element={<Dashboard />} />
				<Route path="/lab" element={<Lab />} />
				<Route path="/library" element={<Library />} />
				<Route path="/profile" element={<Profile />} />
			</Routes>
		</AppShell>
	);
}
