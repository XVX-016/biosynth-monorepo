import React from 'react';
import { Link, useLocation } from 'react-router-dom';

type NavbarProps = {
	onToggleMenu?: () => void;
};

export default function Navbar({ onToggleMenu }: NavbarProps) {
	const location = useLocation();
	const isActive = (to: string, exact = false) =>
		exact ? location.pathname === to : location.pathname.startsWith(to);
	// Workaround for type resolution issues in monorepo
	// eslint-disable-next-line @typescript-eslint/no-explicit-any
	const RLink = (props: any) => React.createElement(Link as any, props);
	return (
		<header className="bg-panel border-b border-aluminum-DEFAULT">
			<div className="mx-auto max-w-7xl px-4 sm:px-6 lg:px-8">
				<div className="h-16 flex items-center justify-between">
					<div className="flex items-center gap-3">
						<button
							type="button"
							className="sm:hidden inline-flex items-center justify-center rounded-lg p-2 text-text-secondary hover:text-text-primary hover:bg-aluminum-light focus:outline-none focus:ring-2 focus:ring-accent-blue"
							aria-label="Toggle navigation"
							onClick={onToggleMenu}
						>
							<svg width="20" height="20" viewBox="0 0 24 24" fill="none" aria-hidden="true">
								<path d="M4 6h16M4 12h16M4 18h16" stroke="currentColor" strokeWidth="2" strokeLinecap="round" />
							</svg>
						</button>
						<span className="text-xl font-bold text-text-primary">BioSynth AI</span>
					</div>
					<nav className="hidden sm:flex items-center gap-2">
						<RLink
							to="/"
							className={`px-3 py-2 rounded-lg font-medium transition-colors ${
								isActive('/', true) ? 'bg-accent-blue text-white' : 'text-text-secondary hover:text-text-primary'
							}`}
						>
							Dashboard
						</RLink>
						<RLink
							to="/lab"
							className={`px-3 py-2 rounded-lg font-medium transition-colors ${
								isActive('/lab') ? 'bg-accent-blue text-white' : 'text-text-secondary hover:text-text-primary'
							}`}
						>
							Lab
						</RLink>
						<RLink
							to="/library"
							className={`px-3 py-2 rounded-lg font-medium transition-colors ${
								isActive('/library') ? 'bg-accent-blue text-white' : 'text-text-secondary hover:text-text-primary'
							}`}
						>
							Library
						</RLink>
						<RLink
							to="/profile"
							className={`px-3 py-2 rounded-lg font-medium transition-colors ${
								isActive('/profile') ? 'bg-accent-blue text-white' : 'text-text-secondary hover:text-text-primary'
							}`}
						>
							Profile
						</RLink>
					</nav>
					<div className="flex items-center">
						<div className="w-8 h-8 rounded-full bg-aluminum-dark/50 border border-aluminum-DEFAULT" aria-label="Profile avatar" />
					</div>
				</div>
			</div>
		</header>
	);
}


