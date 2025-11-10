import React, { useEffect, useState } from 'react';
import { listMolecules } from '../lib/api';

export default function Profile() {
	const [count, setCount] = useState(0);
	useEffect(() => {
		let cancelled = false;
		(async () => {
			try {
				const res = await listMolecules({ limit: 100 });
				if (!cancelled) setCount(res.items?.length ?? 0);
			} catch {
				// ignore for now
			}
		})();
		return () => {
			cancelled = true;
		};
	}, []);

	return (
		<div className="grid grid-cols-12 gap-4">
			<div className="col-span-12 lg:col-span-4">
				<div className="bg-panel rounded-xl shadow-soft border border-aluminum-DEFAULT p-6">
					<div className="flex items-center gap-4">
						<div className="w-16 h-16 rounded-full bg-aluminum-dark/50 border border-aluminum-DEFAULT" />
						<div>
							<h2 className="text-xl font-semibold">Your Profile</h2>
							<p className="text-text-secondary">Local profile</p>
						</div>
					</div>
					<div className="mt-6">
						<label className="block text-sm text-text-secondary mb-1">Display Name</label>
						<input
							type="text"
							placeholder="Your name"
							className="w-full rounded-lg border border-aluminum-DEFAULT bg-aluminum-light px-3 py-2 outline-none focus:ring-2 focus:ring-accent-blue"
						/>
					</div>
				</div>
			</div>
			<div className="col-span-12 lg:col-span-8">
				<div className="bg-panel rounded-xl shadow-soft border border-aluminum-DEFAULT p-6">
					<h3 className="text-lg font-semibold">Summary</h3>
					<div className="mt-4 grid grid-cols-2 gap-4">
						<div className="rounded-xl border border-aluminum-DEFAULT bg-aluminum-light p-4">
							<div className="text-sm text-text-secondary">Saved molecules</div>
							<div className="text-2xl font-bold">{count}</div>
						</div>
						<div className="rounded-xl border border-aluminum-DEFAULT bg-aluminum-light p-4">
							<div className="text-sm text-text-secondary">Recent activity</div>
							<div className="text-2xl font-bold">{Math.min(count, 5)}</div>
						</div>
					</div>
				</div>
			</div>
		</div>
	);
}


