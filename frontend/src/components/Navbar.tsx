import React, { useState } from 'react';
import { Link, useLocation, useNavigate } from 'react-router-dom';
import { useAuthStore } from '../store/authStore';
import Button from './ui/Button';

type NavbarProps = {
	onToggleMenu?: () => void;
};

export default function Navbar({ onToggleMenu }: NavbarProps) {
	const location = useLocation();
	const navigate = useNavigate();
	const { user, signOut, openAuthModal } = useAuthStore();
	const [dropdownOpen, setDropdownOpen] = useState(false);
	const isActive = (to: string, exact = false) =>
		exact ? location.pathname === to : location.pathname.startsWith(to);
	// Workaround for type resolution issues in monorepo
	// eslint-disable-next-line @typescript-eslint/no-explicit-any
	const RLink = (props: any) => React.createElement(Link as any, props);

	const handleSignOut = async () => {
		await signOut();
		setDropdownOpen(false);
		// Don't redirect - user can stay on current page
	};

	const handleSignIn = () => {
		openAuthModal('signin');
	};

	return (
		<header className="bg-white border-b border-lightGrey shadow-neon">
			<div className="mx-auto max-w-7xl px-4 sm:px-6 lg:px-8">
				<div className="h-16 flex items-center justify-between">
					{/* Logo */}
					<div className="flex items-center gap-3">
						<span className="text-xl font-bold text-black">MolForge</span>
					</div>
					
					{/* Center Navigation */}
					<nav className="hidden sm:flex items-center gap-1">
						<RLink
							to="/"
							className={`px-4 py-2 rounded-lg font-medium transition-all duration-200 ${
								isActive('/', true) 
									? 'bg-black text-white' 
									: 'text-darkGrey hover:text-black hover:bg-offwhite'
							}`}
						>
							Home
						</RLink>
						<RLink
							to="/library"
							className={`px-4 py-2 rounded-lg font-medium transition-all duration-200 ${
								isActive('/library') 
									? 'bg-black text-white' 
									: 'text-darkGrey hover:text-black hover:bg-offwhite'
							}`}
						>
							Library
						</RLink>
						<RLink
							to="/lab"
							className={`px-4 py-2 rounded-lg font-medium transition-all duration-200 ${
								isActive('/lab') 
									? 'bg-black text-white' 
									: 'text-darkGrey hover:text-black hover:bg-offwhite'
							}`}
						>
							Lab
						</RLink>
						<RLink
							to="/models"
							className={`px-4 py-2 rounded-lg font-medium transition-all duration-200 ${
								isActive('/models') 
									? 'bg-black text-white' 
									: 'text-darkGrey hover:text-black hover:bg-offwhite'
							}`}
						>
							Models
						</RLink>
						<RLink
							to="/phase10"
							className={`px-4 py-2 rounded-lg font-medium transition-all duration-200 ${
								isActive('/phase10') 
									? 'bg-black text-white' 
									: 'text-darkGrey hover:text-black hover:bg-offwhite'
							}`}
						>
							Phase 10
						</RLink>
						<RLink
							to="/docs"
							className={`px-4 py-2 rounded-lg font-medium transition-all duration-200 ${
								isActive('/docs') 
									? 'bg-black text-white' 
									: 'text-darkGrey hover:text-black hover:bg-offwhite'
							}`}
						>
							Docs
						</RLink>
					</nav>
					
					{/* Right side: Search + Auth */}
					<div className="flex items-center gap-3">
						{/* Search input - hidden on mobile */}
						<div className="hidden md:block">
							<input
								type="text"
								placeholder="Search..."
								className="w-48 rounded-lg border border-lightGrey bg-white text-black px-3 py-1.5 text-sm outline-none focus:ring-2 focus:ring-darkGrey/20 focus:border-darkGrey placeholder:text-midGrey"
							/>
						</div>
						{/* Mobile menu button */}
						<button
							type="button"
							className="sm:hidden inline-flex items-center justify-center rounded-lg p-2 text-darkGrey hover:text-black hover:bg-offwhite focus:outline-none focus:ring-2 focus:ring-darkGrey/20"
							aria-label="Toggle navigation"
							onClick={onToggleMenu}
						>
							<svg width="20" height="20" viewBox="0 0 24 24" fill="none" aria-hidden="true">
								<path d="M4 6h16M4 12h16M4 18h16" stroke="currentColor" strokeWidth="2" strokeLinecap="round" />
							</svg>
						</button>
						{/* Auth buttons */}
						{user ? (
							<div className="relative">
								<button
									onClick={() => setDropdownOpen(!dropdownOpen)}
									className="w-8 h-8 rounded-full bg-lightGrey border border-lightGrey hover:border-darkGrey transition-colors flex items-center justify-center"
									aria-label="Profile menu"
									title={user.email || 'Profile'}
								>
									<svg width="16" height="16" viewBox="0 0 24 24" fill="none" className="text-darkGrey">
										<path d="M12 12c2.21 0 4-1.79 4-4s-1.79-4-4-4-4 1.79-4 4 1.79 4 4 4zm0 2c-2.67 0-8 1.34-8 4v2h16v-2c0-2.66-5.33-4-8-4z" fill="currentColor" />
									</svg>
								</button>
								{dropdownOpen && (
									<>
										<div
											className="fixed inset-0 z-40"
											onClick={() => setDropdownOpen(false)}
										/>
										<div className="absolute right-0 mt-2 w-48 bg-white rounded-lg shadow-lg border border-lightGrey z-50 py-1">
											<Link
												to="/profile"
												onClick={() => setDropdownOpen(false)}
												className="block px-4 py-2 text-sm text-darkGrey hover:bg-offwhite"
											>
												Profile
											</Link>
											<Link
												to="/my-molecules"
												onClick={() => setDropdownOpen(false)}
												className="block px-4 py-2 text-sm text-darkGrey hover:bg-offwhite"
											>
												My Molecules
											</Link>
											<button
												onClick={handleSignOut}
												className="block w-full text-left px-4 py-2 text-sm text-darkGrey hover:bg-offwhite"
											>
												Sign out
											</button>
										</div>
									</>
								)}
							</div>
						) : (
							<Button
								variant="primary"
								size="sm"
								onClick={handleSignIn}
								className="hidden sm:inline-flex"
							>
								Sign in
							</Button>
						)}
					</div>
				</div>
			</div>
			
			{/* Mobile menu */}
			{onToggleMenu && (
				<div className="sm:hidden border-t border-lightGrey bg-white">
					<nav className="px-4 py-2 space-y-1">
						<RLink
							to="/"
							onClick={onToggleMenu}
							className={`block px-3 py-2 rounded-lg transition-all duration-200 ${
								isActive('/', true) 
									? 'bg-black text-white' 
									: 'text-darkGrey hover:text-black hover:bg-offwhite'
							}`}
						>
							Home
						</RLink>
						<RLink
							to="/library"
							onClick={onToggleMenu}
							className={`block px-3 py-2 rounded-lg transition-all duration-200 ${
								isActive('/library') 
									? 'bg-black text-white' 
									: 'text-darkGrey hover:text-black hover:bg-offwhite'
							}`}
						>
							Library
						</RLink>
						<RLink
							to="/lab"
							onClick={onToggleMenu}
							className={`block px-3 py-2 rounded-lg transition-all duration-200 ${
								isActive('/lab') 
									? 'bg-black text-white' 
									: 'text-darkGrey hover:text-black hover:bg-offwhite'
							}`}
						>
							Lab
						</RLink>
						<RLink
							to="/models"
							onClick={onToggleMenu}
							className={`block px-3 py-2 rounded-lg transition-all duration-200 ${
								isActive('/models') 
									? 'bg-black text-white' 
									: 'text-darkGrey hover:text-black hover:bg-offwhite'
							}`}
						>
							Models
						</RLink>
						<RLink
							to="/docs"
							onClick={onToggleMenu}
							className={`block px-3 py-2 rounded-lg transition-all duration-200 ${
								isActive('/docs') 
									? 'bg-black text-white' 
									: 'text-darkGrey hover:text-black hover:bg-offwhite'
							}`}
						>
							Docs
						</RLink>
					</nav>
				</div>
			)}
		</header>
	);
}
