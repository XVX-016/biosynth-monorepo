import React, { useState } from 'react';
import Navbar from '../components/Navbar';
import { Link, useLocation } from 'react-router-dom';

type AppShellProps = {
	children: React.ReactNode;
};

export default function AppShell({ children }: AppShellProps) {
	const [mobileOpen, setMobileOpen] = useState(false);
	const location = useLocation();
	const isActive = (to: string, exact = false) =>
		exact ? location.pathname === to : location.pathname.startsWith(to);
	// Workaround for type resolution issues in monorepo
	// eslint-disable-next-line @typescript-eslint/no-explicit-any
	const RLink = (props: any) => React.createElement(Link as any, props);

	return (
		<div className="min-h-screen bg-ionBlack text-ivory">
			<Navbar onToggleMenu={() => setMobileOpen((v) => !v)} />
			<div className="flex">
				<aside
					className={`${
						mobileOpen ? 'block' : 'hidden'
					} sm:block w-64 shrink-0`}
					aria-label="Sidebar"
				>
					<div className="sticky top-0 sm:top-16 sm:h-[calc(100vh-4rem)] frosted-glass border-r border-chrome/20 shadow-glass backdrop-blur-md">
						<nav className="p-4 space-y-1">
							<RLink
								to="/"
								onClick={() => setMobileOpen(false)}
								className={`block px-3 py-2 rounded-lg transition-all duration-200 ${
									isActive('/', true) 
										? 'bg-plasma-neon text-ionBlack shadow-neon-sm' 
										: 'text-chrome hover:text-ivory hover:bg-frostedGlass'
								}`}
							>
								Dashboard
							</RLink>
							<RLink
								to="/lab"
								onClick={() => setMobileOpen(false)}
								className={`block px-3 py-2 rounded-lg transition-all duration-200 ${
									isActive('/lab') 
										? 'bg-plasma-neon text-ionBlack shadow-neon-sm' 
										: 'text-chrome hover:text-ivory hover:bg-frostedGlass'
								}`}
							>
								Lab
							</RLink>
							<RLink
								to="/library"
								onClick={() => setMobileOpen(false)}
								className={`block px-3 py-2 rounded-lg transition-all duration-200 ${
									isActive('/library') 
										? 'bg-plasma-neon text-ionBlack shadow-neon-sm' 
										: 'text-chrome hover:text-ivory hover:bg-frostedGlass'
								}`}
							>
								Library
							</RLink>
							<RLink
								to="/profile"
								onClick={() => setMobileOpen(false)}
								className={`block px-3 py-2 rounded-lg transition-all duration-200 ${
									isActive('/profile') 
										? 'bg-plasma-neon text-ionBlack shadow-neon-sm' 
										: 'text-chrome hover:text-ivory hover:bg-frostedGlass'
								}`}
							>
								Profile
							</RLink>
						</nav>
					</div>
				</aside>
				<main className="flex-1 p-4 sm:p-6 lg:p-8">{children}</main>
			</div>
		</div>
	);
}


