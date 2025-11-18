import React, { useState } from 'react';
import Navbar from '../components/Navbar';

type AppShellProps = {
	children: React.ReactNode;
};

export default function AppShell({ children }: AppShellProps) {
	const [mobileOpen, setMobileOpen] = useState(false);

	return (
		<div className="min-h-screen bg-white text-black">
			<Navbar onToggleMenu={() => setMobileOpen((v) => !v)} />
			<main className="flex-1 p-4 sm:p-6 lg:p-8">{children}</main>
		</div>
	);
}
